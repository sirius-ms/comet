/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schilller University.
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

package de.unijena.bioinf.projectspace;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.FormulaScore;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.chem.RetentionTime;
import de.unijena.bioinf.ChemistryBase.ms.DetectedAdducts;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.*;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.stream.Collectors;

public final class CompoundContainerId extends ProjectSpaceContainerId {
    public enum Flag {COMPUTING}
    public static final String RANKING_KEY = "rankingScoreType";

    //transient fields
    protected final transient ReentrantReadWriteLock containerLock = new ReentrantReadWriteLock(); //lock for IO
    protected final transient Lock flagsLock = new ReentrantLock(); //lock for non persistent changes
    protected final transient EnumSet<Flag> flags = EnumSet.noneOf(Flag.class);

    // ID defining fields
    private final int compoundIndex;
    @NotNull
    private String directoryName;
    @NotNull
    private String compoundName;

    /**
     * feature id (might be external) of this aligned feature
     */
    @Nullable
    private String featureId;

    /**
     * ID of the compound this feature belongs to
     * ATTENTION: this means compound in the sense of a group of adduct ions that come from the same compound.
     * The term CompoundId refers to an aligned feature -,-
     */
    @Nullable
    private String groupId;
    @Nullable
    private RetentionTime groupRt;

    // fields for fast compound filtering
    @Nullable
    private Double ionMass = null;
    @Nullable
    private PrecursorIonType ionType = null;
    @Nullable
    private DetectedAdducts possibleAdducts = null;
    @NotNull
    private List<Class<? extends FormulaScore>> rankingScores = Collections.emptyList();
    @Nullable
    private RetentionTime rt;
    @Nullable
    private Double confidenceScore;

    protected CompoundContainerId(@NotNull String directoryName, @NotNull String compoundName, int compoundIndex) {
        this(directoryName, compoundName, compoundIndex, null, null, null, null, null, null, null);
    }

    protected CompoundContainerId(@NotNull String directoryName, @NotNull String compoundName, int compoundIndex,
                                  @Nullable Double ionMass, @Nullable PrecursorIonType ionType,
                                  @Nullable RetentionTime rt, @Nullable Double confidenceScore,
                                  @Nullable String featureId, @Nullable String groupId,  @Nullable RetentionTime groupRt) {
        this.directoryName = directoryName;
        this.compoundName = compoundName;
        this.compoundIndex = compoundIndex;
        this.ionMass = ionMass;
        this.ionType = ionType;
        this.rt = rt;
        this.confidenceScore = confidenceScore;
        this.featureId = featureId;
        this.groupId = groupId;
    }

    public boolean hasFlag(Flag flag) {
        return flags.contains(flag);
    }

    public int getCompoundIndex() {
        return compoundIndex;
    }

    @NotNull
    public String getDirectoryName() {
        return directoryName;
    }

    @NotNull
    public String getCompoundName() {
        return compoundName;
    }

    @NotNull
    public Optional<Double> getIonMass() {
        return Optional.ofNullable(ionMass);
    }

    public void setIonMass(@Nullable Double ionMass) {
        this.ionMass = ionMass;
    }

    @NotNull
    public Optional<PrecursorIonType> getIonType() {
        return Optional.ofNullable(ionType);
    }

    public void setIonType(@Nullable PrecursorIonType ionType) {
        this.ionType = ionType;
    }

    public List<Class<? extends FormulaScore>> getRankingScoreTypes() {
        return rankingScores;
    }

    @NotNull
    public Optional<DetectedAdducts> getDetectedAdducts() {
        return Optional.ofNullable(possibleAdducts);
    }

    public void setDetectedAdducts(@Nullable DetectedAdducts possibleAdducts) {
        this.possibleAdducts = possibleAdducts;
    }

    @NotNull
    public Optional<RetentionTime> getRt() {
        return Optional.ofNullable(rt);
    }

    public void setRt(@Nullable RetentionTime rt) {
        this.rt = rt;
    }

    public Optional<Double> getConfidenceScore() {
        return Optional.ofNullable(confidenceScore);
    }

    public void setConfidenceScore(@Nullable Double confidenceScore) {
        this.confidenceScore = confidenceScore;
    }

    @NotNull
    public Optional<String> getFeatureId() {
        return Optional.ofNullable(featureId);
    }

    public void setFeatureId(@Nullable String featureId) {
        this.featureId = featureId;
    }

    public Optional<String> getGroupId() {
        return Optional.ofNullable(groupId);
    }

    public void setGroupId(@Nullable String groupId) {
        this.groupId = groupId;
    }

    @NotNull
    public Optional<RetentionTime> getGroupRt() {
        return Optional.ofNullable(groupRt);
    }

    public void setGroupRt(@Nullable RetentionTime rt) {
        this.groupRt = rt;
    }

    @SafeVarargs
    public final void setRankingScoreTypes(@NotNull Class<? extends FormulaScore>... rankingScores) {
        setRankingScoreTypes(Arrays.asList(rankingScores));
    }

    public void setRankingScoreTypes(@NotNull List<Class<? extends FormulaScore>> rankingScores) {
        this.rankingScores = new ArrayList<>(rankingScores);
    }

    /**
     * This operation is only allowed to be called with careful synchronization within the project space
     */
    void rename(String newName, String newDirName) {
        this.compoundName = newName;
        this.directoryName = newDirName;
    }

    @Override
    public String toString() {
        return directoryName + "@" + getIonMass().map(Math::round).map(String::valueOf).orElse("N/A") + "m/z";
    }

    public Map<String, String> asKeyValuePairs() {
        Map<String, String> kv = new LinkedHashMap<>(3);
        kv.put("index", String.valueOf(getCompoundIndex()));
        kv.put("name", getCompoundName());
        getIonMass().ifPresent(im -> kv.put("ionMass", String.valueOf(im)));
        getIonType().ifPresent(it -> kv.put("ionType", it.toString()));
        getDetectedAdducts().ifPresent(pa -> kv.put("detectedAdducts", pa.toString()));
        getRt().ifPresent(rt -> kv.put("rt", RetentionTime.asStringValue(rt)));
        getConfidenceScore().ifPresent(cs -> kv.put("confidenceScore", String.valueOf(cs)));
        getFeatureId().ifPresent(fid -> kv.put("featureId", featureId));
        getGroupId().ifPresent(fid -> kv.put("groupId", groupId));
        getGroupRt().ifPresent(rt -> kv.put("groupRt", RetentionTime.asStringValue(rt)));


        if (!rankingScores.isEmpty())
            kv.put(RANKING_KEY, rankingScores.stream().map(Score::simplify).collect(Collectors.joining(",")));

        return kv;
    }

    public void setAllNonFinal(final CompoundContainerId cid) {
        if (cid == null || cid == this)
            return;
        setRankingScoreTypes(cid.rankingScores);
        setIonMass(cid.ionMass);
        setIonType(cid.ionType);
        setDetectedAdducts(cid.possibleAdducts);
        setRt(cid.rt);
        setConfidenceScore(cid.confidenceScore);
        setFeatureId(cid.featureId);
        setGroupId(cid.groupId);
        setGroupRt(cid.groupRt);
    }
}
