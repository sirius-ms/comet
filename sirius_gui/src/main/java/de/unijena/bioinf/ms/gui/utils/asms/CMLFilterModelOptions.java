package de.unijena.bioinf.ms.gui.utils.asms;

import lombok.Getter;

@Getter
public class CMLFilterModelOptions {

    private static final CMLFilterModelOptions DISABLED = new CMLFilterModelOptions(null, null, null, 2, 5, 2, 10, 10, false);

    //Parameters for the filtering of measured combinatorial molecule libraries (Affinity selection MS):
    // io
    private final String pathToBBFile;
    private final String scaffoldMf;
    private final String matchedPeaksOutputFilePath;

    // ms1 filter and peak matching filter options:
    private final int minMatchingPeaks;
    private final int numTopPeaks;
    private final double ms1Deviation;
    private final double ms2Deviation;
    private final int numAllowedHydrogenShifts;
    private final boolean isPeakMatchingFilterEnabled;


    public CMLFilterModelOptions(String pathToBBFile, String scaffoldMf, String matchedPeaksOutputFilePath, int minMatchingPeaks, int numTopPeaks, int numAllowedHydrogenShifts, double ms1Deviation, double ms2Deviation, boolean isPeakMatchingFilterEnabled){
        this.pathToBBFile = pathToBBFile;
        this.scaffoldMf = scaffoldMf;
        this.matchedPeaksOutputFilePath = matchedPeaksOutputFilePath;
        this.minMatchingPeaks = minMatchingPeaks;
        this.numTopPeaks = numTopPeaks;
        this.numAllowedHydrogenShifts = numAllowedHydrogenShifts;
        this.ms1Deviation = ms1Deviation;
        this.ms2Deviation = ms2Deviation;
        this.isPeakMatchingFilterEnabled = isPeakMatchingFilterEnabled;
    }

    public static CMLFilterModelOptions disabled(){
        return DISABLED;
    }

    public boolean isMs1FilterActive(){
        return this.pathToBBFile != null && !this.pathToBBFile.isEmpty();
    }

    public boolean isMs2FilterActive(){
        return this.isPeakMatchingFilterEnabled();
    }
}
