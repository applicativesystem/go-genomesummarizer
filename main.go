package main

/*

Author Gaurav Sablok
Date: 2024-11-14
Universitat Potsdam

A golang implementation for summarizing the genome as well as metagenome annotations.
It prepares the ready to use information including all the stats implemented on the
genome as well as metagenome annotations.


*/

import (
	"bufio"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

func main() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatal(err)
	}
	os.Exit(1)
}

var genomesummarize string

var rootCmd = &cobra.Command{
	Use:  "genomesummarize",
	Long: "Run the genome analyzer function and summarize the genome",
	Run:  genomeanalyzeFunc,
}

func init() {
	rootCmd.Flags().
		StringVarP(&genomesummarize, "genomesummarize", "G", "genomeannotate", "summarize the genome")
}

func genomeanalyzeFunc(cmd *cobra.Command, args []string) {
	mRNAID := []string{}

	f, err := os.Open(genomesummarize)
	if err != nil {
		log.Fatal(err)
	}
	fread := bufio.NewScanner(f)
	for fread.Scan() {
		line := fread.Text()
		if strings.Split(string(line), " ")[2] == "mRNA" {
			mRNAprocess := strings.Split(string(line), " ")[8]
			mRNAadapt := strings.Split(mRNAprocess, ";")[0]
			mRNAcapture := strings.Split(mRNAadapt, "=")[2]
			mRNAID = append(mRNAID, mRNAcapture)
		}
	}

	uniqueID := []string{}
	for i := 0; i <= len(mRNAID)-1; i++ {
		if mRNAID[i] == mRNAID[i+1] {
			continue
		} else {
			uniqueID = append(uniqueID, mRNAID[i])
		}
	}

	type informationStruct struct {
		AnnotateType string
		Start        []string
		End          []string
		Estimate     []int
	}

	proteinCap := []informationStruct{}
	cdsCap := []informationStruct{}
	exonCap := []informationStruct{}
	threeUTRCap := []informationStruct{}
	fiveUTRCap := []informationStruct{}
	// mRNACap := []informationStruct{}

	opengff, err := os.Open(genomesummarize)
	if err != nil {
		log.Fatal(err)
	}

	opengffRead := bufio.NewScanner(opengff)

	for opengffRead.Scan() {
		line := opengffRead.Text()
		for j := range uniqueID {
			if strings.Split(string(line), " ")[2] == "CDS" &&
				strings.Split(strings.Split(string(line), " ")[8], ",")[0] == uniqueID[j] {
				cdsAppend := informationStruct{
					AnnotateType: uniqueID[j],
				}
				cdsAppend.Start = append(
					cdsAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				cdsAppend.End = append(
					cdsAppend.End,
					strings.Split(string(line), " ")[4],
				)
				convertStart, _ := strconv.Atoi(strings.Split(string(line), " ")[3])
				convertEnd, _ := strconv.Atoi(strings.Split(string(line), " ")[4])
				cdsAppend.Estimate = append(
					cdsAppend.Estimate, convertEnd-convertStart,
				)
				cdsCap = append(cdsCap, cdsAppend)
			}
		}
	}

	for opengffRead.Scan() {
		line := opengffRead.Text()
		for j := range uniqueID {
			if strings.Split(string(line), " ")[2] == "exon" &&
				strings.Split(string(line), " ")[8] == uniqueID[j] {
				exonAppend := informationStruct{
					AnnotateType: uniqueID[j],
				}
				exonAppend.Start = append(
					exonAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				exonAppend.End = append(
					exonAppend.End, strings.Split(string(line), " ")[4],
				)
				convertStart, _ := strconv.Atoi(strings.Split(string(line), " ")[3])
				convertEnd, _ := strconv.Atoi(strings.Split(string(line), " ")[4])
				exonAppend.Estimate = append(
					exonAppend.Estimate, convertEnd-convertStart,
				)
				exonCap = append(exonCap, exonAppend)
			}
		}
	}

	for opengffRead.Scan() {
		line := opengffRead.Text()
		for j := range uniqueID {
			capturefirst := strings.Split(strings.Split(string(line), " ")[8], ";")[0]
			capturesecond := strings.Split(string(capturefirst), "=")[1]
			capturethird := strings.Split(string(capturesecond), "-")[0]
			if strings.Split(string(line), " ")[2] == "protein" && capturethird == uniqueID[j] {
				proteinAppend := informationStruct{
					AnnotateType: uniqueID[j],
				}
				proteinAppend.Start = append(
					proteinAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				proteinAppend.End = append(
					proteinAppend.End,
					strings.Split(string(line), " ")[4],
				)
				convertStart, _ := strconv.Atoi(strings.Split(string(line), " ")[3])
				convertEnd, _ := strconv.Atoi(strings.Split(string(line), " ")[4])
				proteinAppend.Estimate = append(
					proteinAppend.Estimate, convertEnd-convertStart,
				)
				proteinCap = append(proteinCap, proteinAppend)
			}
		}
	}

	for opengffRead.Scan() {
		line := opengffRead.Text()
		for j := range uniqueID {
			if strings.Split(string(line), " ")[2] == "five_prime_UTR" &&
				strings.Split(string(line), " ")[8] == uniqueID[j] {
				fiveAppend := informationStruct{
					AnnotateType: uniqueID[j],
				}
				fiveAppend.Start = append(
					fiveAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				fiveAppend.End = append(
					fiveAppend.End, strings.Split(string(line), " ")[4],
				)
				convertStart, _ := strconv.Atoi(strings.Split(string(line), " ")[3])
				convertEnd, _ := strconv.Atoi(strings.Split(string(line), " ")[4])
				fiveAppend.Estimate = append(
					fiveAppend.Estimate, convertEnd-convertStart,
				)
				fiveUTRCap = append(fiveUTRCap, fiveAppend)
			}
		}
	}

	for opengffRead.Scan() {
		line := opengffRead.Text()
		for j := range uniqueID {
			if strings.Split(string(line), " ")[2] == "three_prime_UTR" &&
				strings.Split(string(line), " ")[8] == uniqueID[j] {
				threeAppend := informationStruct{
					AnnotateType: uniqueID[j],
				}
				threeAppend.Start = append(
					threeAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				threeAppend.End = append(
					threeAppend.End, strings.Split(string(line), " ")[4],
				)
				convertStart, _ := strconv.Atoi(strings.Split(string(line), " ")[3])
				convertEnd, _ := strconv.Atoi(strings.Split(string(line), " ")[4])
				threeAppend.Estimate = append(
					threeAppend.Estimate, convertEnd-convertStart,
				)
				threeUTRCap = append(threeUTRCap, threeAppend)
			}
		}
	}

	type cummulativeInfo struct {
		id          string
		cummulative int
		length      int
		mean        float64
	}

	proteinInfo := []cummulativeInfo{}
	cdsInfo := []cummulativeInfo{}
	threeInfo := []cummulativeInfo{}
	fiveInfo := []cummulativeInfo{}
	exonInfo := []cummulativeInfo{}

	for i := range proteinCap {
		proteinInfo = append(proteinInfo, cummulativeInfo{
			id:          proteinCap[i].AnnotateType,
			cummulative: sum(proteinCap[i].Estimate),
			length:      len(proteinCap[i].Estimate),
			mean: float64(
				float64(sum(proteinCap[i].Estimate)) / float64(len(proteinCap[i].Estimate)),
			),
		})
	}

	for i := range cdsCap {
		cdsInfo = append(cdsInfo, cummulativeInfo{
			id:          cdsCap[i].AnnotateType,
			cummulative: sum(cdsCap[i].Estimate),
			length:      len(cdsCap[i].Estimate),
			mean: float64(
				float64(sum(cdsCap[i].Estimate)) / float64(len(cdsCap[i].Estimate)),
			),
		})
	}

	for i := range threeUTRCap {
		threeInfo = append(threeInfo, cummulativeInfo{
			id:          threeUTRCap[i].AnnotateType,
			cummulative: sum(threeUTRCap[i].Estimate),
			length:      len(threeUTRCap[i].Estimate),
			mean: float64(
				float64(sum(threeUTRCap[i].Estimate)) / float64(len(threeUTRCap[i].Estimate)),
			),
		})
	}

	for i := range fiveUTRCap {
		fiveInfo = append(fiveInfo, cummulativeInfo{
			id:          fiveUTRCap[i].AnnotateType,
			cummulative: sum(fiveUTRCap[i].Estimate),
			length:      len(fiveUTRCap[i].Estimate),
			mean: float64(
				float64(sum(fiveUTRCap[i].Estimate)) / float64(len(fiveUTRCap[i].Estimate)),
			),
		})
	}

	for i := range exonCap {
		exonInfo = append(exonInfo, cummulativeInfo{
			id:          exonCap[i].AnnotateType,
			cummulative: sum(exonCap[i].Estimate),
			length:      len(exonCap[i].Estimate),
			mean: float64(
				float64(sum(exonCap[i].Estimate)) / float64(len(exonCap[i].Estimate)),
			),
		})
	}
}

func sum(input []int) int {
	sumstart := int(0)
	for i := range input {
		sumstart += input[i]
	}
	return sumstart
}
