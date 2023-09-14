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

package de.unijena.bioinf.chemdb.nitrite;

import de.unijena.bioinf.ChemistryBase.fp.CdkFingerprintVersion;
import de.unijena.bioinf.ChemistryBase.fp.FingerprintVersion;
import de.unijena.bioinf.chemdb.ChemicalNoSQLDatabase;
import de.unijena.bioinf.storage.db.nosql.nitrite.NitriteDatabase;
import org.dizitart.no2.Document;

import java.io.IOException;
import java.nio.file.Path;

public class ChemicalNitriteDatabase extends ChemicalNoSQLDatabase<Document> {

    private final String name;

    public ChemicalNitriteDatabase(Path file) throws IOException {
        super(new NitriteDatabase(file, initMetadata(CdkFingerprintVersion.getDefault())));
        this.name = file.getFileName().toString();
    }

    public ChemicalNitriteDatabase(Path file, FingerprintVersion version) throws IOException {
        super(new NitriteDatabase(file, initMetadata(version)));
        this.name = file.getFileName().toString();
    }

    public NitriteDatabase getStorage(){
        return (NitriteDatabase) storage;
    }

    @Override
    public String getName() {
        return name;
    }
}