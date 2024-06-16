package de.unijena.bioinf.ms.gui.utils.asms;

import lombok.Getter;

import java.util.Optional;

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


    private String currentPathToBBFile;
    private String currentScaffoldMf;
    private String currentMatchedPeaksOutputFilePath;

    private int currentMinMatchingPeaks;
    private int currentNumTopPeaks;
    private int currentNumAllowedHydrogenShifts;
    private double currentMs2Deviation;
    private double currentMs1Deviation;
    private boolean currentIsPeakMatchingFilterEnabled;


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

        this.currentPathToBBFile = pathToBBFile;
        this.currentScaffoldMf = scaffoldMf;
        this.currentMatchedPeaksOutputFilePath = matchedPeaksOutputFilePath;
        this.currentMinMatchingPeaks = minMatchingPeaks;
        this.currentNumTopPeaks = numTopPeaks;
        this.currentNumAllowedHydrogenShifts = numAllowedHydrogenShifts;
        this.currentMs1Deviation = ms1Deviation;
        this.currentMs2Deviation = ms2Deviation;
        this.currentIsPeakMatchingFilterEnabled = isPeakMatchingFilterEnabled;
    }

    public static CMLFilterModelOptions disabled(){
        return DISABLED;
    }

    public boolean isMs1FilterActive(){
        return this.currentPathToBBFile != null && !Optional.ofNullable(this.pathToBBFile).orElse("").equals(this.currentPathToBBFile);
    }

    public boolean isMs2FilterActive(){
        return this.isCurrentIsPeakMatchingFilterEnabled();
    }

}
