/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer and Sebastian Böcker,
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
 *  You should have received a copy of the GNU Lesser General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.ms.persistence.model.core.feature;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import com.fasterxml.jackson.annotation.JsonIgnore;
import de.unijena.bioinf.ms.persistence.model.core.run.LCMSRun;
import de.unijena.bioinf.ms.persistence.model.core.spectrum.IsotopePattern;
import de.unijena.bioinf.ms.persistence.model.core.spectrum.MSData;
import jakarta.persistence.Id;
import lombok.*;
import lombok.experimental.SuperBuilder;

import javax.annotation.Nullable;
import java.util.List;
import java.util.Optional;

import static com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility.ANY;
import static com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility.NONE;

/**
 * Features aligned over several {@link LCMSRun}s (same m/z and RT)
 * Annotations/identifications are calculated and annotated at this level
 */
@Getter
@Setter
@NoArgsConstructor
@AllArgsConstructor
@SuperBuilder
@ToString(callSuper = true)
@JsonAutoDetect(fieldVisibility = ANY, getterVisibility = NONE, setterVisibility = NONE, isGetterVisibility = NONE)
public class AlignedFeatures extends AbstractAlignedFeatures {

    @Id
    private long alignedFeatureId;

    /**
     * ID of the compound this aligned feature belongs to
     */
    private long compoundId;

    /**
     * Detected adducts with score and source information for this feature
     */
    private DetectedAdducts detectedAdducts;

    @JsonIgnore
    @ToString.Exclude
    private List<AlignedIsotopicFeatures> isotopicFeatures;

    public Optional<List<AlignedIsotopicFeatures>> getIsotopicFeatures() {
        return Optional.ofNullable(isotopicFeatures);
    }


    public static AlignedFeatures singleton(Feature feature) {
        return singleton(feature, (MSData) null);
    }
    public static AlignedFeatures singleton(Feature feature, @Nullable IsotopePattern isotopePattern) {
        // TODO add also MS/MS spectra (and merged spectrum?)
        MSData msData1 = isotopePattern != null ? MSData.builder().isotopePattern(isotopePattern).build() : null;
        return singleton(feature, msData1);
    }
    public static AlignedFeatures singleton(Feature feature, @Nullable MSData msData) {
        return AlignedFeatures.builder()
                .features(List.of(feature))
                .averageMass(feature.averageMass)
                .apexMass(feature.apexMass)
                .apexIntensity(feature.apexIntensity)
                .snr(feature.snr)
                .retentionTime(feature.retentionTime)
                .ionType(feature.getIonType())
                .msData(msData)
                .build();
    }


}
