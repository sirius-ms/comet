package de.unijena.bioinf.io;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.datastructures.*;

import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;

public class BuildingBlockReader {

    public static BuildingBlock[][] readBuildingBlocks(File file) throws IOException, UnknownElementException {
        try(BufferedReader reader = Files.newBufferedReader(file.toPath())){
            return readBuildingBlocks(reader);
        }
    }

    public static BuildingBlock[][] readBuildingBlocks(Reader r) throws IOException, UnknownElementException {
        try(BufferedReader reader = new BufferedReader(r)){
            return readBuildingBlocks(reader);
        }
    }

    public static BuildingBlock[][] readBuildingBlocks(BufferedReader reader) throws IOException, UnknownElementException {
        reader.readLine(); // this is just the header of the CSV file
        final ArrayList<ArrayList<BuildingBlock>> buildingBlocks = new ArrayList<>();
        String currentLine = reader.readLine();

        while(currentLine != null){
            final String[] row = currentLine.split(",");
            final int currentPos = Integer.parseInt(row[0]);
            if(currentPos >= buildingBlocks.size()){
                final int diff = currentPos - buildingBlocks.size() + 1;
                for(int i = 0; i < diff; i++) buildingBlocks.add(new ArrayList<>());
            }

            final BuildingBlock bb = new BuildingBlock(getFinalMolecularFormula(row[2], row[3]), row[1]);
            final ArrayList<BuildingBlock> bbs = buildingBlocks.get(currentPos);
            final int bbIdx = Integer.parseInt(row[4])-1; // the indices in the file start with 1
            if(bbIdx >= bbs.size()){
                final int diff = bbIdx - bbs.size() + 1;
                for(int i = 0; i < diff; i++) bbs.add(null);
            }
            bbs.set(bbIdx, bb);
            currentLine = reader.readLine();
        }

        // Transform to BuildingBlock matrix:
        return buildingBlocks.stream().map(bbList -> bbList.toArray(BuildingBlock[]::new)).toArray(BuildingBlock[][]::new);
    }

    private static MolecularFormula getFinalMolecularFormula(String mf, String loss) throws UnknownElementException {
        MolecularFormula molFormula = MolecularFormula.parse(mf);
        MolecularFormula lossMf = MolecularFormula.parse(loss);
        return molFormula.subtract(lossMf);
    }

}
