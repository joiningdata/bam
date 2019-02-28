package bam

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"
)

// An Index contains information to allow fast lookup
// of sequences aligning to a region of reference sequence.
type Index struct {
	Refs []IndexReference
}

// A IndexReference contains alignment info for the reference sequence.
type IndexReference struct {
	// Bins group aligned sequences into a tree structure.
	Bins map[uint32]Bin

	// Intervals have the linear index of aligned sequences.
	Intervals []Offset

	// Unmapped reads are placed into a single Chunk.
	Unmapped Chunk

	// TotalMapped read-segments for this reference.
	TotalMapped uint64

	// TotalUnmapped read-segments for this reference.
	TotalUnmapped uint64
}

// A Bin contains a list of Chunks.
type Bin []Chunk

// Chunk is a part of a Bin, representing a section of alignments.
type Chunk struct{ Begin, End Offset }

// Offset represents a virtual offset within the compressed BAM file.
type Offset uint64

// Compressed returns the offset of the compressed block start.
func (o Offset) Compressed() int64 {
	return int64(o >> 16)
}

// Uncompressed returns the offset within the uncompressed block.
func (o Offset) Uncompressed() uint16 {
	return uint16(o & 0xFFFF)
}

/////////

// LoadIndex for a BAM file. You probably don't need this, it will automatically be loaded
// via the Load("file.bam") method as long as "file.bam.bai" exists.
func LoadIndex(filename string) (*Index, error) {
	le := binary.LittleEndian

	ff, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	f := &Index{}
	tmp := make([]byte, 8)
	_, err = io.ReadFull(ff, tmp)
	if err != nil {
		return nil, err
	}
	if bytes.Equal(tmp, []byte{'B', 'A', 'M', 1}) {
		return nil, fmt.Errorf("bam: invalid index file '%v'", tmp)
	}
	n := int32(le.Uint32(tmp[4:]))
	f.Refs = make([]IndexReference, n)
	for i, r := range f.Refs {
		_, err = io.ReadFull(ff, tmp[:4])
		if err != nil {
			return nil, err
		}
		nb := int32(le.Uint32(tmp[:4]))
		r.Bins = make(map[uint32]Bin, nb)

		BAMProgressFunc(float64(i*100) / float64(n))

		for j := int32(0); j < nb; j++ {
			_, err = io.ReadFull(ff, tmp)
			if err != nil {
				return nil, err
			}
			bid := le.Uint32(tmp[:4])
			nc := int32(le.Uint32(tmp[4:]))
			b := make([]Chunk, nc)
			err = binary.Read(ff, le, &b)
			if err != nil {
				return nil, err
			}
			if bid == 37450 {
				// Unmapped reads are held/recorded separately
				r.Unmapped = b[0]
				r.TotalMapped = uint64(b[1].Begin)
				r.TotalUnmapped = uint64(b[1].End)
				continue
			}
			r.Bins[bid] = b
		}

		_, err = io.ReadFull(ff, tmp[:4])
		if err != nil {
			return nil, err
		}
		ni := int32(le.Uint32(tmp[:4]))
		r.Intervals = make([]Offset, ni)
		err = binary.Read(ff, le, &r.Intervals)
		if err != nil {
			return nil, err
		}
		f.Refs[i] = r
	}
	BAMProgressFunc(-1.0)
	return f, err
}

func (r *IndexReference) getBin(beginPos, endPos uint64) uint32 {
	endPos = (endPos - 1) >> 14
	beginPos >>= 14

	if beginPos == endPos {
		return ((1<<15)-1)/7 + uint32(beginPos)
	}
	if (beginPos >> 3) == (endPos >> 3) {
		return ((1<<12)-1)/7 + uint32(beginPos>>3)
	}
	if (beginPos >> 6) == (endPos >> 6) {
		return ((1<<9)-1)/7 + uint32(beginPos>>6)
	}
	if (beginPos >> 9) == (endPos >> 9) {
		return ((1<<6)-1)/7 + uint32(beginPos>>9)
	}
	if (beginPos >> 12) == (endPos >> 12) {
		return ((1<<3)-1)/7 + uint32(beginPos>>12)
	}
	return 0
}

func (r *IndexReference) getBins(beginPos, endPos uint64) []uint32 {
	res := make([]uint32, 1, ((1<<18)-1)/7)

	endPos = (endPos - 1) >> 14
	beginPos >>= 14

	for k := 1 + beginPos>>12; k <= 1+(endPos>>12); k++ {
		res = append(res, uint32(k))
	}
	for k := 9 + beginPos>>9; k <= 1+(endPos>>9); k++ {
		res = append(res, uint32(k))
	}
	for k := 73 + beginPos>>6; k <= 73+(endPos>>6); k++ {
		res = append(res, uint32(k))
	}
	for k := 585 + beginPos>>3; k <= 585+(endPos>>3); k++ {
		res = append(res, uint32(k))
	}
	for k := 4681 + beginPos; k <= 4681+endPos; k++ {
		res = append(res, uint32(k))
	}
	return res
}
