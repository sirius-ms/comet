/*
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2015 Kai Dührkop
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.unijena.bioinf.fteval;

import net.sf.jniinchi.INCHI_RET;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.similarity.Tanimoto;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class ChemicalSimilarity {

    private final List<Compound> compounds;
    private final List<IFingerprinter> fingerprinters;
    private final EvalDB db;

    public ChemicalSimilarity(EvalDB db) {
        this(db, new PubchemFingerprinter(DefaultChemObjectBuilder.getInstance()), new MACCSFingerprinter());
    }

    public ChemicalSimilarity(EvalDB db, IFingerprinter... fingerprinter) {
        this.compounds = new ArrayList<Compound>();
        this.fingerprinters = new ArrayList<IFingerprinter>(Arrays.asList(fingerprinter));
        this.db = db;

    }

    public void addInchi(String name, String inchi) {
        final BitSet[] bitsets = new BitSet[fingerprinters.size()];
        try {
            InChIToStructure str = InChIGeneratorFactory.getInstance().getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
            if (str.getReturnStatus() != INCHI_RET.OKAY) {
                System.err.println("Error for InChI:\n" + inchi);
                System.err.println(str.getMessage());
                return;
            }
            IAtomContainer mol = str.getAtomContainer();

            int k = 0;
            for (IFingerprinter finger : fingerprinters) {
                bitsets[k++] = finger.getBitFingerprint(mol).asBitSet();
            }
            compounds.add(new Compound(db.removeExtName(new File(name)), bitsets, inchi));
        } catch (CDKException e) {
            throw new RuntimeException(e);
        }
    }

    public void add(String name) throws IOException {
        InChIGeneratorFactory factory = null;
        try {
            factory = InChIGeneratorFactory.getInstance();
        } catch (CDKException e) {
            throw new RuntimeException(e);
        }
        final File sdf = db.sdf(name);
        final BufferedReader reader = new BufferedReader(new FileReader(sdf));
        final ISimpleChemObjectReader chemReader = new ReaderFactory().createReader(reader);
        try {
            final AtomContainer mol = chemReader.read(new AtomContainer());
            final BitSet[] bitsets = new BitSet[fingerprinters.size()];
            int k = 0;
            for (IFingerprinter finger : fingerprinters) {
                bitsets[k++] = finger.getBitFingerprint(mol).asBitSet();
            }
            final InChIGenerator gen = factory.getInChIGenerator(mol);
            if (gen.getReturnStatus() != INCHI_RET.OKAY) {
                System.err.println(gen.getMessage());
            }
            compounds.add(new Compound(db.removeExtName(new File(name)), bitsets, gen.getInchi()));
        } catch (CDKException e) {
            throw new RuntimeException(e);
        }
    }

    public List<Compound> getCompounds() {
        return compounds;
    }

    public float[][][] computeTanimoto() {
        final String[] names = new String[compounds.size()];
        int k = 0;
        Collections.sort(compounds, new Comparator<Compound>() {
            @Override
            public int compare(Compound o1, Compound o2) {
                return o1.name.compareTo(o2.name);
            }
        });
        for (Compound c : compounds) names[k++] = c.name;
        final float[][][] matrix = new float[fingerprinters.size()][compounds.size()][compounds.size()];
        for (int f = 0; f < fingerprinters.size(); ++f) {
            for (int row = 0; row < names.length; ++row) {
                for (int col = row + 1; col < names.length; ++col) {
                    try {
                        matrix[f][row][col] = matrix[f][col][row] = Tanimoto.calculate(compounds.get(row).fingerprints[f], compounds.get(col).fingerprints[f]);
                    } catch (CDKException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
        }
        return matrix;
    }

    public File[] getDirectories() {
        final File[] names = new File[fingerprinters.size()];
        int k = 0;
        for (IFingerprinter finger : fingerprinters) {
            final Class<? extends IFingerprinter> klass = finger.getClass();
            final String name = klass.getSimpleName().replace("Fingerprinter", "");
            final File dir = db.fingerprint(name).getParentFile();
            names[k++] = dir;
        }
        return names;
    }

    public static class Compound {
        String inchi;
        String name;
        private BitSet[] fingerprints;

        private Compound(String name, BitSet[] fingerprints, String inchi) {
            this.name = name;
            this.fingerprints = fingerprints;
            this.inchi = inchi;
        }
    }
}
