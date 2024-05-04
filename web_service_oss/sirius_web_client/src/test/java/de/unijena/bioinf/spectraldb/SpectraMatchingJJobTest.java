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

package de.unijena.bioinf.spectraldb;

import com.github.f4b6a3.tsid.Tsid;
import de.unijena.bioinf.ChemistryBase.jobs.SiriusJobs;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.utils.OrderedSpectrumDelegate;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums;
import de.unijena.bioinf.chemdb.ChemicalDatabaseException;
import de.unijena.bioinf.spectraldb.entities.Ms2ReferenceSpectrum;
import de.unijena.bionf.spectral_alignment.CosineQuerySpectrum;
import de.unijena.bionf.spectral_alignment.CosineQueryUtils;
import de.unijena.bionf.spectral_alignment.SpectralMatchingType;
import de.unijena.bionf.spectral_alignment.SpectralMatchMasterJJob;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Rule;
import org.junit.Test;
import org.mockito.junit.MockitoJUnit;
import org.mockito.junit.MockitoRule;
import org.mockito.quality.Strictness;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class SpectraMatchingJJobTest {

    @Rule public MockitoRule mockito = MockitoJUnit.rule().strictness(Strictness.STRICT_STUBS);

    @Test
    public void testGetCosineQueries() {
        Deviation dev = new Deviation(10);
        SpectraMatchingJJob job = new SpectraMatchingJJob(null, null);
        CosineQueryUtils utils = new CosineQueryUtils(SpectralMatchingType.MODIFIED_COSINE.getScorer(new Deviation(10)));

        Ms2Spectrum<Peak> q1 = new MutableMs2Spectrum(new SimpleSpectrum(new double[] {1d, 2d}, new double[] {1d, 2d}), 30, null, 0);
        Ms2Spectrum<Peak> q2 = new MutableMs2Spectrum(new SimpleSpectrum(new double[] {1d, 3d}, new double[] {1d, 2d}), 40, null, 0);

        List<CosineQuerySpectrum> cosineQueries = new ArrayList<>();
        cosineQueries.addAll(job.getCosineQueries(utils, dev,  q1.getPrecursorMz(), Arrays.asList(q1,q2)));
        CosineQuerySpectrum cq1 = cosineQueries.get(0);
        CosineQuerySpectrum cq2 = cosineQueries.get(1);

        assertEquals(2, cosineQueries.size());
//        assertTrue(Spectrums.haveEqualPeaks(q1, cq1)); //todo this does not work because spectra are processed and normalized
//        assertTrue(Spectrums.haveEqualPeaks(q2, cq2));
        assertEquals(q1.getPrecursorMz(), cq1.getPrecursorMz(), 1e-9);
        assertEquals(q2.getPrecursorMz(), cq2.getPrecursorMz(), 1e-9);
        assertEquals(0, cq1.getIndex());
        assertEquals(1, cq2.getIndex());
    }

    @Test
    public void testGetJJobsForDb() throws ChemicalDatabaseException {
        Deviation dev = new Deviation(10);
        SpectraMatchingJJob job = new SpectraMatchingJJob(null, null);
        CosineQueryUtils utils = new CosineQueryUtils(SpectralMatchingType.INTENSITY.getScorer(dev));

        Ms2Spectrum<Peak> q1 = new MutableMs2Spectrum();
        Ms2Spectrum<Peak> q2 = new MutableMs2Spectrum();

        CosineQuerySpectrum cq1 = utils.createQueryWithoutLoss(new OrderedSpectrumDelegate<>(q1), 0);
        CosineQuerySpectrum cq2 = utils.createQueryWithoutLoss(new OrderedSpectrumDelegate<>(q2), 0);

        SimpleSpectrum rs1 = SimpleSpectrum.empty();
        SimpleSpectrum rs2 = SimpleSpectrum.empty();

        Ms2ReferenceSpectrum r1 = Ms2ReferenceSpectrum.builder().spectrum(rs1).build();
        Ms2ReferenceSpectrum r2 = Ms2ReferenceSpectrum.builder().spectrum(rs2).build();


        assertNotEquals(q1, q2);
        assertNotEquals(cq1, cq2);
        assertNotEquals(r1, r2);

        final List<Ms2ReferenceSpectrum> references = Arrays.asList(r1, r2);
        assertNotNull(r1.getSpectrum());
        assertNotNull(r2.getSpectrum());

        List<SpectralMatchMasterJJob> jobs = job.getAlignmentJJobs(utils, Arrays.asList(cq1, cq2), references);

        SpectralMatchMasterJJob j1 = jobs.get(0);
        SpectralMatchMasterJJob j2 = jobs.get(1);

        assertEquals(2, jobs.size());
        assertEquals(2, j1.getQueries().size());
        assertEquals(2, j2.getQueries().size());

        Pair<CosineQuerySpectrum, CosineQuerySpectrum> p11 = j1.getQueries().get(0);
        Pair<CosineQuerySpectrum, CosineQuerySpectrum> p12 = j1.getQueries().get(1);
        Pair<CosineQuerySpectrum, CosineQuerySpectrum> p21 = j2.getQueries().get(0);
        Pair<CosineQuerySpectrum, CosineQuerySpectrum> p22 = j2.getQueries().get(1);

        assertEquals(cq1, p11.getLeft());
        assertEquals(cq1, p12.getLeft());
        assertEquals(cq2, p21.getLeft());
        assertEquals(cq2, p22.getLeft());

        assertNull(r1.getSpectrum());
        assertNull(r2.getSpectrum());
    }

    @Test
    public void testExtractResults() {
        SpectraMatchingJJob job = new SpectraMatchingJJob(null, null);
        CosineQueryUtils utils = new CosineQueryUtils(SpectralMatchingType.INTENSITY.getScorer(new Deviation(10)));

        SimpleSpectrum query = new SimpleSpectrum(new double[] {1d}, new double[] {1d});
        int queryIndex = 5;
        CosineQuerySpectrum cosineQuery = utils.createQueryWithoutLoss(new IndexedQuerySpectrumWrapper(query, queryIndex), 1);

        Ms2ReferenceSpectrum emptyReference = Ms2ReferenceSpectrum.builder()
                .spectrum(SimpleSpectrum.empty())
                .precursorMz(0d)
                .build();
        CosineQuerySpectrum cosineEmptyReference = utils.createQueryWithoutLoss(new Ms2ReferenceSpectrumWrapper(emptyReference), 0);

        String libraryId = "library id";
        String libraryName = "library name";
        long UUID = Tsid.fast().toLong();

        SimpleSpectrum referenceSpectrum = new SimpleSpectrum(new double[] {1d, 2d}, new double[] {0.5d, 1d});
        Ms2ReferenceSpectrum reference = Ms2ReferenceSpectrum.builder()
                .spectrum(referenceSpectrum)
                .precursorMz(3d)
                .libraryId(libraryId)
                .libraryName(libraryName)
                .uuid(UUID)
                .build();
        List<Ms2ReferenceSpectrum> refs = List.of(reference);
        CosineQuerySpectrum cosineReference = utils.createQueryWithoutLoss(new Ms2ReferenceSpectrumWrapper(reference), 3);
        cosineReference.setIndex(0);

        SpectralMatchMasterJJob matchJob = new SpectralMatchMasterJJob(utils, Arrays.asList(Pair.of(cosineQuery, cosineEmptyReference), Pair.of(cosineQuery, cosineReference)));
        matchJob.setClearInput(false);
        SiriusJobs.getGlobalJobManager().submitJob(matchJob);

        List<SpectralSearchResult.SearchResult> results = job.extractResults(matchJob, refs).toList();

        assertEquals(1, results.size());
        SpectralSearchResult.SearchResult result = results.get(0);

        assertEquals(queryIndex, result.getQuerySpectrumIndex());
        assertEquals(libraryId, result.getDbId());
        assertEquals(libraryName, result.getDbName());
        assertEquals(UUID, result.getUuid());
        assertEquals(utils.cosineProduct(cosineQuery, cosineReference), result.getSimilarity());
    }
}