package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/joiningdata/bam"
)

func commas(n int) string {
	s := fmt.Sprint(n)
	r := len(s) % 3
	if r == 0 {
		r = 3
	}
	cs := []string{s[:r]}
	for j := r; j < len(s); j += 3 {
		cs = append(cs, s[j:j+3])
	}
	return strings.Join(cs, ",")
}

func main() {
	maxmem := flag.String("m", "500M", "maximum memory size to use")
	listRefs := flag.Bool("l", false, "list reference sequence info")
	listBins := flag.Bool("lb", false, "list bin details for each reference sequence")
	refName := flag.String("r", "", "query named reference only")
	startPos := flag.Int64("s", 0, "start position for alignment map (0-based)")
	endPos := flag.Int64("e", -1, "end position for alignment map (0-based)")
	flag.Parse()

	*maxmem = strings.ToUpper(*maxmem)
	*maxmem = strings.TrimSuffix(*maxmem, "B")
	if *maxmem != "500M" {
		lc := (*maxmem)[len(*maxmem)-1]
		var mult int64
		switch lc {
		case 'K':
			mult = 1024
		case 'M':
			mult = 1024 * 1024
		case 'G':
			mult = 1024 * 1024 * 1024
		case 'T':
			mult = 1024 * 1024 * 1024 * 1024
		default:
			mult = 1
		}
		if mult > 1 {
			*maxmem = (*maxmem)[:len(*maxmem)-1]
		}
		num, err := strconv.ParseInt(*maxmem, 10, 64)
		if err != nil {
			log.Println(err)
		} else {
			bam.MaxBAMMemory = num * mult
		}
	}

	b, err := bam.Load(flag.Arg(0))
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(1)
	}

	refID := -1
	if *refName != "" {
		for i, r := range b.References {
			if *refName == r.Name {
				refID = i

				if *endPos < 0 {
					*endPos = int64(r.Length)
				}

				break
			}
		}
	}

	if *listRefs || *listBins || refID == -1 {
		for i, r := range b.References {
			if i != refID && refID != -1 {
				continue
			}
			fmt.Printf("%4d: %-30s   %12s bp\n", i+1, r.Name, commas(r.Length))

			if b.Index == nil || !*listBins {
				continue
			}
			bins := b.Index.Refs[i]
			for j, b := range bins.Bins {
				if len(b) == 1 {
					fmt.Printf("      Bin %5d: %d-%d\n", j, b[0].Begin, b[0].End)
				} else if len(b) == 2 {
					fmt.Printf("      Bin %5d: %d-%d and %d-%d\n", j, b[0].Begin, b[0].End, b[1].Begin, b[1].End)
				} else {
					fmt.Printf("      Bin %5d: %d-%d and %d-%d + %d more\n", j, b[0].Begin, b[0].End, b[1].Begin, b[1].End, len(b)-2)
				}
			}
		}
	}

	if refID == -1 {
		os.Exit(0)
	}

	data := b.GetMap(int32(refID), uint64(*startPos), uint64(*endPos))
	for _, row := range data {
		fmt.Println(len(row), row)
	}
}
