/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer and Sebastian Böcker,
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
 *  You should have received a copy of the GNU Lesser General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.ms.rest.model;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonInclude;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusJobInput;
import de.unijena.bioinf.ms.rest.model.covtree.CovtreeJobInput;
import de.unijena.bioinf.ms.rest.model.fingerid.FingerprintJobInput;
import de.unijena.bioinf.ms.rest.model.msnovelist.MsNovelistJobInput;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

@JsonInclude(JsonInclude.Include.NON_NULL)
@JsonIgnoreProperties(ignoreUnknown = true)
public class JobInputs {
    private List<FingerprintJobInput> fingerprintJobInputs = new ArrayList<>();
    private List<CanopusJobInput> canopusJobInputs = new ArrayList<>();
    private List<CovtreeJobInput> covtreeJobInputs = new ArrayList<>();
    private List<MsNovelistJobInput> msnovelistJobInputs = new ArrayList<>();

    public List<FingerprintJobInput> getFingerprintJobInputs() {
        return Collections.unmodifiableList(fingerprintJobInputs);
    }

    @JsonIgnore
    public boolean hasFingerprintJobs(){
        return fingerprintJobInputs != null && !fingerprintJobInputs.isEmpty();
    }

    public void setFingerprintJobInputs(List<FingerprintJobInput> fingerprintJobInputs) {
        this.fingerprintJobInputs = fingerprintJobInputs;
    }


    @JsonIgnore
    public void addFingerprintJobInput(FingerprintJobInput fingerprintJobInput) {
        addFingerprintJobInputs(List.of(fingerprintJobInput));
    }

    @JsonIgnore
    public void addFingerprintJobInputs(List<FingerprintJobInput> fingerprintJobInputs) {
        if (this.fingerprintJobInputs == null)
            this.fingerprintJobInputs = new ArrayList<>();
        this.fingerprintJobInputs.addAll(fingerprintJobInputs);
    }

    public List<CanopusJobInput> getCanopusJobInputs() {
        return Collections.unmodifiableList(canopusJobInputs);
    }

    @JsonIgnore
    public boolean hasCanopusJobs(){
        return canopusJobInputs != null && !canopusJobInputs.isEmpty();
    }

    public void setCanopusJobInputs(List<CanopusJobInput> canopusJobInputs) {
        this.canopusJobInputs = canopusJobInputs;
    }

    @JsonIgnore
    public void addCanopusJobInput(CanopusJobInput canopusJobInputs) {
        addCanopusJobInputs(List.of(canopusJobInputs));
    }

    @JsonIgnore
    public void addCanopusJobInputs(List<CanopusJobInput> canopusJobInputs) {
        if (this.canopusJobInputs == null)
            this.canopusJobInputs = new ArrayList<>();
        this.canopusJobInputs.addAll(canopusJobInputs);
    }

    public List<CovtreeJobInput> getCovtreeJobInputs() {
        return Collections.unmodifiableList(covtreeJobInputs);
    }

    @JsonIgnore
    public boolean hasCovtreeJobs(){
        return covtreeJobInputs != null && !covtreeJobInputs.isEmpty();
    }

    public void setCovtreeJobInputs(List<CovtreeJobInput> covtreeJobInputs) {
        this.covtreeJobInputs = covtreeJobInputs;
    }

    @JsonIgnore
    public void addCovtreeJobInput(CovtreeJobInput covtreeJobInputs) {
        addCovtreeJobInputs(List.of(covtreeJobInputs));
    }

    @JsonIgnore
    public void addCovtreeJobInputs(List<CovtreeJobInput> covtreeJobInputs) {
        if (this.covtreeJobInputs == null)
            this.covtreeJobInputs = new ArrayList<>();
        this.covtreeJobInputs.addAll(covtreeJobInputs);
    }

    public List<MsNovelistJobInput> getMsNovelistJobInputs() {
        return Collections.unmodifiableList(msnovelistJobInputs);
    }

    @JsonIgnore
    public boolean hasMsNovelistJobs(){
        return msnovelistJobInputs != null && !msnovelistJobInputs.isEmpty();
    }

    public void setMsNovelistSSJobInputs(List<MsNovelistJobInput> msnovelistJobInputs) {
        this.msnovelistJobInputs = msnovelistJobInputs;
    }

    @JsonIgnore
    public void addMsNovelistJobInput(MsNovelistJobInput msnovelistJobInputs) {
        addMsNovelistJobInputs(List.of(msnovelistJobInputs));
    }

    @JsonIgnore
    public void addMsNovelistJobInputs(List<MsNovelistJobInput> msnovelistJobInputs) {
        if (this.msnovelistJobInputs == null)
            this.msnovelistJobInputs = new ArrayList<>();
        this.msnovelistJobInputs.addAll(msnovelistJobInputs);
    }

    @JsonIgnore
    public Map<JobTable, List<?>> asMap() {
        return Map.of(
                JobTable.JOBS_FINGERID, fingerprintJobInputs,
                JobTable.JOBS_CANOPUS, canopusJobInputs,
                JobTable.JOBS_COVTREE, covtreeJobInputs,
                JobTable.JOBS_MSNOVELIST, msnovelistJobInputs
        );
    }
    @JsonIgnore
    public void addJobInput(Object jobInput, JobTable type) {
        switch (type) {
            case JOBS_FINGERID -> addFingerprintJobInput(((FingerprintJobInput) jobInput));
            case JOBS_CANOPUS -> addCanopusJobInput((CanopusJobInput) jobInput);
            case JOBS_COVTREE -> addCovtreeJobInput((CovtreeJobInput) jobInput);
            case JOBS_MSNOVELIST -> addMsNovelistJobInput((MsNovelistJobInput) jobInput);
        }
    }

    @JsonIgnore
    public boolean hasJobs(){
        return hasCanopusJobs() || hasFingerprintJobs() || hasCovtreeJobs() || hasMsNovelistJobs();
    }

    @JsonIgnore
    public int size(){
        return fingerprintJobInputs.size() + canopusJobInputs.size() + covtreeJobInputs.size() + msnovelistJobInputs.size();
    }
}
