package bam

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"strconv"
	"strings"
)

var bgzfEOF = []byte{
	0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff,
	6, 0, 0x42, 0x43, 2, 0, 0x1b, 0, 3, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
}

// An AlignmentMap represents a sequence alignment/map.
type AlignmentMap struct {
	blockCache map[int64][]byte

	Index *Index

	Header     string
	References []Reference

	Alignments []*bamAlignment
}

type bgzfByteReader struct {
	*os.File
	x [1]byte
}

func (f *bgzfByteReader) ReadByte() (byte, error) {
	_, err := f.Read(f.x[:])
	return f.x[0], err
}

// Load a BAM dataset from the file.
func Load(filename string) (*AlignmentMap, error) {
	ff, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	f := &AlignmentMap{}
	f.blockCache = make(map[int64][]byte)

	/////////
	// check for proper End-of-file marker
	_, err = ff.Seek(-int64(len(bgzfEOF)), io.SeekEnd)
	if err != nil {
		return nil, err
	}
	tmp := make([]byte, len(bgzfEOF))
	_, err = io.ReadFull(ff, tmp)
	if err != nil {
		return nil, err
	}
	if !bytes.Equal(tmp, bgzfEOF) {
		return nil, fmt.Errorf("invalid end-of-file marker (possibly truncated?)")
	}
	ff.Seek(0, io.SeekStart)
	/////////
	bff := bufio.NewReader(ff)
	zr, err := gzip.NewReader(bff)
	if err != nil {
		ff.Close()
		return nil, err
	}

	var remainder []byte

	truepos := int64(0)
	for {
		zr.Multistream(false)
		h := zr.Header
		if h.Extra[0] != 'B' || h.Extra[1] != 'C' {
			panic("not a BAM file (invalid subfield id)")
		}
		if 2 != binary.LittleEndian.Uint16(h.Extra[2:]) {
			panic("not a BAM file (invalid subfield length)")
		}
		bsize := binary.LittleEndian.Uint16(h.Extra[4:]) + 1

		/// read the data here
		data, err := ioutil.ReadAll(zr)
		if err != nil {
			zr.Close()
			ff.Close()
			return nil, err
		}

		f.blockCache[truepos] = data

		if len(f.blockCache) == 1 {
			// parse the header + initial block
			remainder = f.parseHead(data[:])
		} else {
			// copy the partial block from the last chunk to
			// the beginning of this one.
			newchunk := make([]byte, len(remainder))
			copy(newchunk, remainder)
			newchunk = append(newchunk, data...)
			remainder = f.parseNext(newchunk)
		}

		// workaround for go bug #30230
		truepos += int64(bsize)
		ff.Seek(truepos, io.SeekStart)
		bff.Reset(ff)

		// move to the next chunk
		err = zr.Reset(bff)
		if err == io.EOF {
			break
		}
	}

	zr.Close()
	ff.Close()
	f.Index, err = LoadIndex(filename + ".bai")
	if os.IsNotExist(err) {
		log.Println("warning: no index available for", filename)
		err = nil
	}
	return f, err
}

// Reference sequence name and length.
type Reference struct {
	// Name of the reference sequence.
	Name string
	// Length of the reference sequence.
	Length int
}

func (b *AlignmentMap) parseHead(r []byte) []byte {
	le := binary.LittleEndian

	headLength := le.Uint32(r[4:])
	b.Header = string(r[8 : 8+headLength])
	numRefs := int(le.Uint32(r[8+headLength:]))

	offs := 12 + int(headLength)
	for i := 0; i < numRefs; i++ {
		br := Reference{}
		nameLength := int(le.Uint32(r[offs:]))
		br.Name = string(r[offs+4 : offs+4+nameLength-1])
		br.Length = int(le.Uint32(r[offs+4+nameLength:]))
		b.References = append(b.References, br)
		offs += 8 + nameLength
	}

	return b.parseNext(r[offs:])
}

func (b *AlignmentMap) parseNext(r []byte) []byte {
	le := binary.LittleEndian

	for len(r) >= 4 {
		blocksize := int(le.Uint32(r))
		if len(r) < (blocksize + 4) {
			break
		}
		ba := parseAlignment(r[4 : 4+blocksize])
		b.Alignments = append(b.Alignments, ba)

		r = r[4+blocksize:]
	}

	return r
}

type bamAlignment struct {
	refID        int32
	pos          int32
	mapq         uint8
	bin          uint16
	cigarOpCount uint16
	flag         uint16
	seqLen       int32
	nextRefID    int32
	nextPos      int32
	tlen         int32

	readName    string
	cigarPacked []uint32
	seqPacked   []uint8
	qual        string

	AuxData map[string]interface{}
}

func parseAlignment(r []byte) *bamAlignment {
	b := &bamAlignment{}
	le := binary.LittleEndian

	b.refID = int32(le.Uint32(r[0:]))
	b.pos = int32(le.Uint32(r[4:]))
	readNameLen := r[8]
	b.mapq = r[9]
	b.bin = le.Uint16(r[10:])
	b.cigarOpCount = le.Uint16(r[12:])
	b.flag = le.Uint16(r[14:])
	b.seqLen = int32(le.Uint32(r[16:]))
	b.nextRefID = int32(le.Uint32(r[20:]))
	b.nextPos = int32(le.Uint32(r[24:]))
	b.tlen = int32(le.Uint32(r[28:]))
	offs := 32 + int(readNameLen)
	b.readName = string(r[32 : offs-1])

	b.cigarPacked = make([]uint32, b.cigarOpCount)
	bb := bytes.NewBuffer(r[offs:])
	binary.Read(bb, le, &b.cigarPacked)
	offs += 4 * int(b.cigarOpCount)
	b.seqPacked = r[offs : offs+(int(1+b.seqLen)/2)]
	if (b.seqLen % 2) == 1 {
		// ensure sequence past end is set to 0
		b.seqPacked[len(b.seqPacked)-1] &= 0xF0
	}
	offs += (int(1+b.seqLen) / 2)
	b.qual = string(r[offs : offs+int(b.seqLen)])
	offs += int(b.seqLen)

	b.AuxData = make(map[string]interface{})
	for offs < len(r) {
		tag := string(r[offs : offs+2])
		vtype := r[offs+2]
		offs += 3
		switch vtype {
		case 'A', 'c', 'C':
			if vtype == 'c' {
				b.AuxData[tag] = int8(r[offs])
			} else {
				b.AuxData[tag] = r[offs]
			}
			offs++
		case 's', 'S':
			x := le.Uint16(r[offs:])
			if vtype == 's' {
				b.AuxData[tag] = int16(x)
			} else {
				b.AuxData[tag] = x
			}
			offs += 2
		case 'i', 'I':
			x := le.Uint32(r[offs:])
			if vtype == 'i' {
				b.AuxData[tag] = int32(x)
			} else {
				b.AuxData[tag] = x
			}
			offs += 4
		case 'f':
			var x float32
			bb = bytes.NewBuffer(r[offs:])
			binary.Read(bb, le, &x)
			b.AuxData[tag] = x
			offs += 4
		case 'Z':
			o := offs
			for r[o] != 0 {
				o++
			}
			b.AuxData[tag] = string(r[offs:o])
			offs = o + 1
		case 'H':
			x := make([]byte, 0, 64)
			o := offs
			for r[o] != 0 {
				z, _ := strconv.ParseUint(string(r[o:o+2]), 16, 8)
				x = append(x, byte(z))
				o += 2
			}
			b.AuxData[tag] = x
			offs = o + 1
		case 'B':
			vtype = r[offs+1]
			count := le.Uint32(r[offs+2:])

			offs += 6
			bb = bytes.NewBuffer(r[offs:])

			var arr interface{}
			switch vtype {
			case 'c':
				arr = make([]int8, count)
				offs += int(count)
			case 'C':
				arr = make([]uint8, count)
				offs += int(count)
			case 's':
				arr = make([]int16, count)
				offs += int(count * 2)
			case 'S':
				arr = make([]uint16, count)
				offs += int(count * 2)
			case 'i':
				arr = make([]int32, count)
				offs += int(count * 4)
			case 'I':
				arr = make([]uint32, count)
				offs += int(count * 4)
			case 'f':
				arr = make([]float32, count)
				offs += int(count * 4)
			}
			binary.Read(bb, le, &arr)
			b.AuxData[tag] = arr

		default:
			log.Printf("aux data type '%c' not implemented", vtype)
		}
	}

	return b
}

// GetMap returns an alignment of the region.
func (b *AlignmentMap) GetMap(refID int32, beginPos, endPos uint64) []string {
	var result []string
	ref := b.References[refID]
	if beginPos > uint64(ref.Length) || endPos > uint64(ref.Length) {
		panic("invalid range")
	}
	iref := b.Index.Refs[refID]
	bid := iref.getBin(beginPos, endPos)
	bin := iref.Bins[bid]

	for _, chunk := range bin {
		p1 := chunk.Begin.Compressed()
		po := chunk.Begin.Uncompressed()
		p2 := chunk.End.Compressed()

		done := false
		var remainder []byte
		for pi := p1; pi <= p2; pi++ {
			r := b.blockCache[pi][po:]
			if len(remainder) > 0 {
				newchunk := make([]byte, len(remainder))
				copy(newchunk, remainder)
				newchunk = append(newchunk, r...)
				r = newchunk
			}
			le := binary.LittleEndian

			for len(r) >= 4 {
				blocksize := int(le.Uint32(r))
				if len(r) < (blocksize + 4) {
					break
				}
				ba := parseAlignment(r[4 : 4+blocksize])
				if ba.refID != refID {
					done = true
					break
				}

				// alignment is actually in range?
				if ba.pos+ba.tlen >= int32(beginPos) &&
					ba.pos <= int32(endPos) {

					seq := UnpackSequence(ba.seqPacked)
					px := int(ba.pos) - int(beginPos)
					pad := ""
					if px > 0 {
						pad = strings.Repeat(" ", px)
					} else {
						px = -px
						if px >= len(seq) {
							seq = ""
						} else {
							seq = seq[px:]
						}
					}
					seq = pad + seq
					epad := int(endPos - beginPos)
					if len(seq) > epad {
						seq = seq[:epad]
					} else {
						seq = seq + strings.Repeat(" ", epad-len(seq))
					}
					result = append(result, seq)
				}
				r = r[4+blocksize:]
			}
			remainder = r
			if done {
				// done with this chunk
				break
			}
			po = 0
		}
	}
	return result
}

func UnpackSequence(packed []byte) string {
	packmap := []byte("=ACMGRSVTWYHKDBN")
	var r []byte
	for _, p := range packed {
		r = append(r, packmap[(p>>4)&0x0F], packmap[p&0x0F])
	}
	// odd number of characters?
	if r[len(r)-1] == '=' {
		r = r[:len(r)-1]
	}
	return string(r)
}