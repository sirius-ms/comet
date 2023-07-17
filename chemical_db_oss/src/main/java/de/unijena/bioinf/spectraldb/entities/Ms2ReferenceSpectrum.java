/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2023 Bright Giant GmbH
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
 *  You should have received a copy of the GNU General Public License along with SIRIUS.
 *  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.spectraldb.entities;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import com.fasterxml.jackson.databind.ser.std.ToStringSerializer;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.CollisionEnergy;
import de.unijena.bioinf.ChemistryBase.ms.MsInstrumentation;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import de.unijena.bioinf.chemdb.DBLink;
import lombok.*;

@Getter
@Setter
@NoArgsConstructor
@AllArgsConstructor
@Builder
public class Ms2ReferenceSpectrum {

    @Builder.Default
    private long id = -1L;

    /**
     * This is the InChiKey (2D) to map spectra to a standardized SIRIUS structure candidate.
     * In rare cases ot might not match the actual smiles of the measured compound due to standardization.
     */
    private String candidateInChiKey;

    /**
     * The adduct / ion type / precursor type
     */
    @JsonSerialize(using = ToStringSerializer.class)
    @JsonDeserialize(using = SimpleSerializers.PrecursorIonTypeDeserializer.class)
    private PrecursorIonType precursorIonType;

    private double precursorMz;

    private double ionMass;

    private int msLevel = 0;

    @JsonSerialize(using = ToStringSerializer.class)
    @JsonDeserialize(using = SimpleSerializers.CollisionEnergyDeserializer.class)
    private CollisionEnergy collisionEnergy;

    @JsonDeserialize(using = SimpleSerializers.MSInstrumentationDeserializer.class)
    private MsInstrumentation instrumentation;

    /**
     * Molecular formula of the measured compound. Must match candidateInChiKey and  smiles
     */
    @JsonSerialize(using = ToStringSerializer.class)
    @JsonDeserialize(using = SimpleSerializers.MolecularFormulaDeserializer.class)
    private MolecularFormula formula;

    private String name;
    /**
     * This is the representation of the structure that produced this spectrum.
     */
    private String smiles;

    private String libraryName;
    /**
     * Identifier of the spectral library, e,g. nist or massbank id
     * Most libraries seem to use splash now anyway.
     */
    private String libraryId;

    private String splash;

    private SimpleSpectrum spectrum;

    @JsonIgnore
    public DBLink getSpectralDbLink() {
        return new DBLink(libraryName, libraryId);
    }


    @JsonIgnore
    public void setSpectralDbLink(DBLink spectralDbLink) {
        this.libraryName = spectralDbLink.name;
        this.libraryId = spectralDbLink.id;
    }
}
