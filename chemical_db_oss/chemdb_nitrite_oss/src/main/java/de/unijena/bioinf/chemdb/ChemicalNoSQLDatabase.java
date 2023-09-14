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

package de.unijena.bioinf.chemdb;

import de.unijena.bioinf.ChemistryBase.chem.InChI;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.fp.FingerprintVersion;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.chemdb.nitrite.serializers.FingerprintCandidateDbSerializer;
import de.unijena.bioinf.chemdb.nitrite.wrappers.FingerprintCandidateWrapper;
import de.unijena.bioinf.chemdb.nitrite.wrappers.FingerprintWrapper;
import de.unijena.bioinf.spectraldb.SpectralNoSQLDatabase;
import de.unijena.bioinf.storage.db.nosql.*;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.MissingResourceException;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import static de.unijena.bioinf.chemdb.ChemDbTags.TAG_DATE;

public abstract class ChemicalNoSQLDatabase<Doctype> extends SpectralNoSQLDatabase<Doctype> implements AbstractChemicalDatabase {

    public static final String SETTINGS_COLLECTION = "CUSTOM-DB-SETTINGS";

    public ChemicalNoSQLDatabase(Database<Doctype> database) throws IOException {
        super(database);
    }

    protected static Metadata initMetadata(FingerprintVersion version) throws IOException {
        Metadata metadata = SpectralNoSQLDatabase.initMetadata();
        return metadata
                // TODO is the FingerprintWrapper actually necessary or does FingerprintCandidateWrapper work for all use cases?
                .addRepository(
                        FingerprintWrapper.class,
                        new Index("inchikey", IndexType.UNIQUE)
                ).addRepository(
                        FingerprintCandidateWrapper.class,
                        "id",
                        new Index("formula", IndexType.NON_UNIQUE),
                        new Index("mass", IndexType.NON_UNIQUE),
                        new Index("candidate.inchikey", IndexType.UNIQUE),
                        new Index("candidate.bitset", IndexType.NON_UNIQUE),
                        new Index("candidate.name", IndexType.NON_UNIQUE)
                ).addCollection(
                        SETTINGS_COLLECTION
                ).setOptionalFields(
                        FingerprintCandidateWrapper.class, "candidate"
                ).addSerialization(
                        FingerprintCandidate.class,
                        new FingerprintCandidateDbSerializer(),
                        new JSONReader.FingerprintCandidateDeserializer(version)
                );
    }

    @Override
    public List<FormulaCandidate> lookupMolecularFormulas(double ionMass, Deviation deviation, PrecursorIonType ionType) throws ChemicalDatabaseException {
        try {
            final double mass = ionType.precursorMassToNeutralMass(ionMass);
            final double from = mass - deviation.absoluteFor(mass);
            final double to = mass + deviation.absoluteFor(mass);
            return this.storage.findStr(Filter.build().and().gte("mass", from).lte("mass", to), FingerprintCandidateWrapper.class, "candidate")
                    .map(c -> c.getCandidate().toFormulaCandidate(ionType)).toList();
        } catch (IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    @Override
    public boolean containsFormula(MolecularFormula formula) throws ChemicalDatabaseException {
        try {
            return this.storage.count(new Filter().eq("formula", formula.toString()), FingerprintCandidateWrapper.class, 0, 1) > 0;
        } catch (IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    @Override
    public List<CompoundCandidate> lookupStructuresByFormula(MolecularFormula formula) throws ChemicalDatabaseException {
        try {
            return storage.findStr(new Filter().eq("formula", formula.toString()), FingerprintCandidateWrapper.class, "candidate")
                    .map(fc -> fc.getCandidate().toCompoundCandidate()).toList();
        } catch (RuntimeException | IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    @Override
    public <T extends Collection<FingerprintCandidate>> T lookupStructuresAndFingerprintsByFormula(MolecularFormula formula, T fingerprintCandidates) throws ChemicalDatabaseException {
        try {
            StreamSupport.stream(
                    storage.find(Filter.build().eq("formula", formula.toString()), FingerprintCandidateWrapper.class, "candidate").spliterator(), false
            ).map(FingerprintCandidateWrapper::getCandidate).forEach(fingerprintCandidates::add);
            return fingerprintCandidates;
        } catch (IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    public Stream<FingerprintCandidate> lookupFingerprintsByInchisStr(Iterable<String> inchi_keys) throws ChemicalDatabaseException {
        try {
            Object[] keys = StreamSupport.stream(inchi_keys.spliterator(), false).toArray();
            return StreamSupport.stream(storage.find(Filter.build().in("candidate.inchikey", keys), FingerprintCandidateWrapper.class, "candidate"
            ).spliterator(), false).map(FingerprintCandidateWrapper::getCandidate);
        } catch (IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    @Override
    public List<FingerprintCandidate> lookupFingerprintsByInchis(Iterable<String> inchi_keys) throws ChemicalDatabaseException {
        return lookupFingerprintsByInchisStr(inchi_keys).toList();
    }

    @Override
    public List<InChI> lookupManyInchisByInchiKeys(Iterable<String> inchi_keys) throws ChemicalDatabaseException {
        return lookupFingerprintsByInchisStr(inchi_keys).map(FingerprintCandidate::getInchi).toList();
    }

    @Override
    public List<FingerprintCandidate> lookupManyFingerprintsByInchis(Iterable<String> inchi_keys) throws ChemicalDatabaseException {
        return this.lookupFingerprintsByInchis(inchi_keys);
    }

    @Override
    public List<FingerprintCandidate> lookupFingerprintsByInchi(Iterable<CompoundCandidate> compounds) throws ChemicalDatabaseException {
        Iterable<String> inchiIt = () -> {
            Iterator<CompoundCandidate> it = compounds.iterator();
            return new Iterator<>() {
                @Override
                public boolean hasNext() {
                    return it.hasNext();
                }

                @Override
                public String next() {
                    return it.next().getInchiKey2D();
                }
            };
        };
        return lookupFingerprintsByInchis(inchiIt);
    }

    @Override
    public List<InChI> findInchiByNames(List<String> names) throws ChemicalDatabaseException {
        try {
            return storage.findStr(new Filter().in("candidate.name", names.toArray()), FingerprintCandidateWrapper.class, "candidate")
                    .map(fc -> fc.getCandidate().getInchi()).toList();
        } catch (IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    @Override
    public String getChemDbDate() throws ChemicalDatabaseException {
        try {
            return storage.findStr(Filter.build().eq("key", TAG_DATE), Tag.class, 0, 1)
                    .map(Tag::getValue).findFirst()
                    .orElseThrow(() -> new MissingResourceException("Could not find Database Date tag", "Tag", TAG_DATE));
        } catch (IOException e) {
            throw new ChemicalDatabaseException(e);
        }
    }

    @Override
    public void annotateCompounds(List<? extends CompoundCandidate> sublist) throws ChemicalDatabaseException {
        throw new UnsupportedOperationException("Compounds of this database are already Annotated! So annotation is not supported!");
    }

    @Override
    public void close() throws IOException {
        this.storage.close();
    }
}
