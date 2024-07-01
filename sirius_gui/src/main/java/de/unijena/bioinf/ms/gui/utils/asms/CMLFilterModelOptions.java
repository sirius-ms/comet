package de.unijena.bioinf.ms.gui.utils.asms;

import lombok.Getter;

import java.util.List;

@Getter
public class CMLFilterModelOptions {

    private static final CMLFilterModelOptions DISABLED = new CMLFilterModelOptions(null, null, null, null,2, 5, 2, 10, 10, false);

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
    private final List<String> fragmentTypes;


    public CMLFilterModelOptions(String pathToBBFile, String scaffoldMf, String matchedPeaksOutputFilePath, List<String> fragmentTypes, int minMatchingPeaks, int numTopPeaks, int numAllowedHydrogenShifts, double ms1Deviation, double ms2Deviation, boolean isPeakMatchingFilterEnabled){
        this.pathToBBFile = pathToBBFile;
        this.scaffoldMf = scaffoldMf;
        this.matchedPeaksOutputFilePath = matchedPeaksOutputFilePath;
        this.fragmentTypes = fragmentTypes;
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

    public boolean isOutputPathActive(){
        return this.matchedPeaksOutputFilePath != null && !this.matchedPeaksOutputFilePath.isEmpty();
    }

    public String getFragmentTypesString(){
        if(this.fragmentTypes == null || this.fragmentTypes.isEmpty()) return null;
        final StringBuilder strBuilder = new StringBuilder();
        strBuilder.append(this.fragmentTypes.get(0));
        for(int i = 1; i < this.fragmentTypes.size(); i++)
            strBuilder.append(",").append(this.fragmentTypes.get(i));
        return strBuilder.toString();
    }
}
