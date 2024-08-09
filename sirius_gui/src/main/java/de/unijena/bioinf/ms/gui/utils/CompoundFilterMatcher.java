package de.unijena.bioinf.ms.gui.utils;/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2021 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schiller University.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

import ca.odell.glazedlists.matchers.Matcher;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.datastructures.BBFragment;
import de.unijena.bioinf.datastructures.CMLMolecule;
import de.unijena.bioinf.sirius.Sirius;
import de.unijena.bioinf.ChemistryBase.chem.FormulaConstraints;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.chem.RetentionTime;
import de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums;
import de.unijena.bioinf.ChemistryBase.ms.utils.WrapperSpectrum;
import de.unijena.bioinf.datastructures.CMLCandidates;
import de.unijena.bioinf.datastructures.OrderedCombinatorialMoleculeLibrary;
import de.unijena.bioinf.ms.gui.properties.GuiProperties;
import de.unijena.bioinf.ms.gui.utils.asms.CMLFilterModelOptions;
import de.unijena.bioinf.ms.nightsky.sdk.model.*;
import de.unijena.bioinf.ms.properties.PropertyManager;
import de.unijena.bioinf.projectspace.FormulaResultBean;
import de.unijena.bioinf.projectspace.InstanceBean;
import de.unijena.bioinf.sirius.ProcessedPeak;
import org.jetbrains.annotations.NotNull;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class CompoundFilterMatcher implements Matcher<InstanceBean> {
    final CompoundFilterModel filterModel;
    private final GuiProperties properties;

    public CompoundFilterMatcher(GuiProperties properties, CompoundFilterModel filterModel) {
        this.filterModel = filterModel;
        this.properties = properties;
    }

    @Override
    public boolean matches(InstanceBean item) {
        double mz = item.getIonMass();
        double rt = item.getRT().map(RetentionTime::getRetentionTimeInSeconds).orElse(Double.NaN);
        //todo hotfix, since the confidence score is a FormulaScore which sets all NaN to -Infinity (after computation, and thus also in project space)
        double confidence = item.getConfidenceScore(properties.getConfidenceDisplayMode()).filter(conf -> !Double.isInfinite(conf)).orElse(Double.NaN);

        {
            if (mz < filterModel.getCurrentMinMz())
                return false;
            if (filterModel.isMaxMzFilterActive() && mz > filterModel.getCurrentMaxMz())
                return false;
        }

        if (!Double.isNaN(rt)) { //never filter NaN because RT is just not available
            if (rt < filterModel.getCurrentMinRt())
                return false;
            if (filterModel.isMaxRtFilterActive() && rt > filterModel.getCurrentMaxRt())
                return false;
        }

        if (!Double.isNaN(confidence)) {
            if (filterModel.isMinConfidenceFilterActive() && confidence < filterModel.getCurrentMinConfidence())
                return false;
            if (filterModel.isMaxConfidenceFilterActive() && confidence > filterModel.getCurrentMaxConfidence())
                return false;
        } else if (filterModel.isMinConfidenceFilterActive()) { // filter NaN if min filter is set
            return false;
        }

        if (filterModel.isHasMs1() && !item.getSourceFeature().isHasMs1())
            return false;

        if (filterModel.isHasMsMs() && !item.getSourceFeature().isHasMsMs())
            return false;

        if (filterModel.isAdductFilterActive() && !filterModel.getAdducts().contains(item.getIonType()))
            return false;

        if (item.getSourceFeature().getQuality() != null) //always allow to pass the filter if now quality data is available
            if (filterModel.getFeatureQualityFilter().isEnabled() && !filterModel.getFeatureQualityFilter().isQualitySelected(item.getSourceFeature().getQuality()))
                return false;

        return anyIOIntenseFilterMatches(item, filterModel);
    }

    private boolean anyIOIntenseFilterMatches(InstanceBean item, CompoundFilterModel filterModel) {
        if (filterModel.getIoQualityFilters().stream().anyMatch(CompoundFilterModel.QualityFilter::isEnabled)) {
            AlignedFeatureQuality qualityReport = item.getQualityReport();
            if (qualityReport != null) { //always allow to pass the filter if now quality data is available
                Map<String, Category> categories = qualityReport.getCategories();
                for (CompoundFilterModel.QualityFilter filter : filterModel.getIoQualityFilters()) {
                    if (filter.isEnabled()) {
                        Category q = categories.get(filter.getName());
                        if (q != null && !filter.isQualitySelected(q.getOverallQuality()))
                            return false;
                    }
                }
            }
        }

        if (filterModel.isElementFilterEnabled())
            if (!matchesElementFilter(item, filterModel)) return false;


        if (filterModel.isLipidFilterEnabled())
            if (!matchesLipidFilter(item, filterModel)) return false;

        if (filterModel.isDbFilterEnabled())
            if (!matchesDBFilter(item, filterModel)) return false;

        if (filterModel.isTagHidingEnabled())
            if(filterModel.featureSubtractionMatches(item)) return false;

        if(filterModel.isCmlMs1FilterActive()){
            if(!matchesCmlMs1Filter(item, filterModel)) return false;

            if(filterModel.isCmlMs2FilterActive())
                if(!matchesCmlMs2Filter(item, filterModel)) return false;
        }

        return true;
    }

    private boolean matchesCmlMs1Filter(InstanceBean item, CompoundFilterModel filterModel) {
        final Ms2Experiment exp = item.asMs2Experiment();
        CMLCandidates candidates = exp.getAnnotationOrNull(CMLCandidates.class);
        if(candidates == null) {
            final OrderedCombinatorialMoleculeLibrary cmlLibrary = filterModel.getCompoundList().getCmlLibrary();
            final CMLFilterModelOptions cmlFilterOptions = filterModel.getCmlFilterOptions();
            candidates = new CMLCandidates(cmlLibrary, exp, cmlFilterOptions.getMs1Deviation());
            exp.setAnnotation(CMLCandidates.class, candidates);
        }

        return !candidates.getCandidates().isEmpty();
    }

    private boolean matchesCmlMs2Filter(InstanceBean item, CompoundFilterModel filterModel) {
        final Ms2Experiment exp = item.asMs2Experiment();
        final PrecursorIonType ionType = item.getIonType().isIonizationUnknown() ? PrecursorIonType.fromString("[M+H]+") : item.getIonType();
        final List<ProcessedPeak> mergedPeaks = new Sirius().preprocessForMs2Analysis(exp).getMergedPeaks();
        mergedPeaks.remove(mergedPeaks.size()-1); // remove the precursor peak
        mergedPeaks.sort(new ProcessedPeak.RelativeIntensityComparator()); // sort according to the relative intensity

        final CMLFilterModelOptions cmlFilterOptions = filterModel.getCmlFilterOptions();
        final List<String> fragmentTypes = cmlFilterOptions.getFragmentTypes();
        final Deviation ms2Deviation = new Deviation(cmlFilterOptions.getMs2Deviation());
        final int minMatchingPeaks = cmlFilterOptions.getMinMatchingPeaks();
        final int numTopPeaks = cmlFilterOptions.getNumTopPeaks();
        final int numHydrogenShifts = cmlFilterOptions.getNumAllowedHydrogenShifts();
        final double HYDROGEN_MASS = MolecularFormula.getHydrogen().getMass();

        final ArrayList<CMLMolecule> candidates = exp.getAnnotation(CMLCandidates.class).orElse(new CMLCandidates(filterModel.getCompoundList().getCmlLibrary(), exp, cmlFilterOptions.getMs1Deviation())).getCandidates();
        final ArrayList<String> selectedOutputMatchedPeaks = new ArrayList<>();
        final int startPeakIdx = mergedPeaks.size() > numTopPeaks ? mergedPeaks.size() - numTopPeaks : 0;
        boolean foundMatchingCompound = false; // true -> there exists at least one candidate whose fragments match at least 'minMatchingPeaks'
        // (at least one candidate structure can explain 'minMatchingPeaks' in the spectrum)

        for(final CMLMolecule candidate : candidates){
            final ArrayList<String> matchedPeaksStrings = new ArrayList<>();
            final List<BBFragment> fragments = fragmentTypes == null ? candidate.createAllBBFragments() : candidate.createSpecificBBFragments(fragmentTypes);
            int matchedPeaks = 0;

            for(int peakIdx = startPeakIdx; peakIdx < mergedPeaks.size(); peakIdx++){
                final ProcessedPeak peak = mergedPeaks.get(peakIdx);
                final double neutralPeakMass = ionType.precursorMassToNeutralMass(peak.getMass());
                boolean isMatched = false;

                for(final BBFragment fragment : fragments){
                    final double fragmentMass = fragment.getMass();
                    for(int h = -numHydrogenShifts; h <= numHydrogenShifts; h++){
                        final double shiftedFragmentMass = fragmentMass + h * HYDROGEN_MASS;
                        if(ms2Deviation.inErrorWindow(shiftedFragmentMass, neutralPeakMass)){
                            if(!isMatched){
                                isMatched = true;
                                matchedPeaks++;
                            }
                            matchedPeaksStrings.add(item.getFeatureId()+","+peak.getMass()+","+peak.getRelativeIntensity()+","+
                                    candidate.getName()+","+fragment.getFragmentTypeString()+","+h);
                        }
                    }
                }
            }

            if(matchedPeaks >= minMatchingPeaks){
                selectedOutputMatchedPeaks.addAll(matchedPeaksStrings);
                foundMatchingCompound = true;
            }
        }

        if(foundMatchingCompound){
            final StringBuilder fileStrBuilder = new StringBuilder();
            final StringBuilder loggerStrBuilder = new StringBuilder();
            for (final String matchedPeaksString : selectedOutputMatchedPeaks) {
                final String[] strArray = matchedPeaksString.split(",");
                loggerStrBuilder.append("[feature_ID: ").append(strArray[0]).append(" | matched peak: (").append(strArray[1])
                        .append(";").append(strArray[2]).append(") | compound: ").append(strArray[3]).append(" | fragment_ID: ")
                        .append(strArray[4]).append(" | H-shifts: ").append(strArray[5]).append("]\n");
                fileStrBuilder.append(matchedPeaksString).append("\n");
            }
            LoggerFactory.getLogger(this.getClass()).info(loggerStrBuilder.toString());

            if(filterModel.isCmlOutputPathActive()) {
                synchronized (this) {
                    try (BufferedWriter fileWriter = new BufferedWriter(new FileWriter(cmlFilterOptions.getMatchedPeaksOutputFilePath(), true))) {
                        fileWriter.write(fileStrBuilder.toString());
                    } catch (IOException e) {
                        LoggerFactory.getLogger(getClass()).warn("Couldn't write information about matched peaks into the file.");
                        e.printStackTrace();
                    }
                }
            }
            return true;
        }else{
            return false;
        }
    }

    private boolean matchesLipidFilter(InstanceBean item, CompoundFilterModel filterModel) {
        boolean hasAnyLipidHit = item.getFormulaCandidates().stream().anyMatch(FormulaResultBean::isLipid);
        return (filterModel.getLipidFilter() == CompoundFilterModel.LipidFilter.ANY_LIPID_CLASS_DETECTED && hasAnyLipidHit)
                || (filterModel.getLipidFilter() == CompoundFilterModel.LipidFilter.NO_LIPID_CLASS_DETECTED && !hasAnyLipidHit);
    }

    private boolean matchesDBFilter(InstanceBean item, CompoundFilterModel filterModel) {
        final int k;
        List<String> filterDbs;
        if (filterModel.isDbFilterEnabled()) {
            k = filterModel.getDbFilter().getNumOfCandidates();
            filterDbs = filterModel.getDbFilter().getDbs().stream().map(SearchableDatabase::getDatabaseId).toList();
        } else {
            k = 1;
            filterDbs = null;
        }

        if (k == 0)
            return false;

        final PageStructureCandidateFormula candidates = item.getStructureCandidatesPage(k, false);

        if (candidates == null || candidates.getContent() == null || candidates.getContent().isEmpty())
            return false;

        if (filterDbs == null)
            return true;

        return candidates.getContent().stream()
                .map(StructureCandidateFormula::getDbLinks)
                .filter(Objects::nonNull).flatMap(List::stream)
                .map(DBLink::getName).distinct()
                .filter(Objects::nonNull)
                .anyMatch(filterDbs::contains);
    }

    private boolean matchesElementFilter(InstanceBean item, CompoundFilterModel filterModel) {
        CompoundFilterModel.ElementFilter filter = filterModel.getElementFilter();
        @NotNull FormulaConstraints constraints = filter.constraints;
        return item.getFormulaAnnotationAsBean().map(fc ->
                (filter.matchFormula && constraints.isSatisfied(fc.getMolecularFormulaObj(), fc.getAdductObj().getIonization()))
                        || (filter.matchPrecursorFormula && constraints.isSatisfied(fc.getAdductObj().neutralMoleculeToMeasuredNeutralMolecule(fc.getMolecularFormulaObj()), fc.getAdductObj().getIonization()))
        ).orElse(false);
    }
}
