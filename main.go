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
	}

	proteinCap := []informationStruct{}
	cdsCap := []informationStruct{}
	exonCap := []informationStruct{}
	threeUTRCap := []informationStruct{}
	fiveUTRCap := []informationStruct{}

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
					AnnotateType: strings.Split(string(line), " ")[2],
				}
				cdsAppend.Start = append(
					cdsAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				cdsAppend.End = append(
					cdsAppend.End,
					strings.Split(string(line), " ")[4],
				)
				cdsCap := append(cdsCap, cdsAppend)
			}
		}
	}

	for opengffRead.Scan() {
		line := opengffRead.Text()
		for j := range uniqueID {
			if strings.Split(string(line), " ")[2] == "exon" &&
				strings.Split(string(line), " ")[8] == uniqueID[j] {
				exonAppend := informationStruct{
					AnnotateType: strings.Split(string(line), " ")[2],
				}
				exonAppend.Start = append(
					exonAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				exonAppend.End = append(
					exonAppend.End, strings.Split(string(line), " ")[4],
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
					AnnotateType: strings.Split(string(line), " ")[2],
				}
				proteinAppend.Start = append(
					proteinAppend.Start,
					strings.Split(string(line), " ")[3],
				)
				proteinAppend.End = append(proteinAppend.End, strings.Split(string(line), " ")[4])
				proteinCap = append(proteinCap, proteinAppend)
			}
		}
	}
}
