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

package de.unijena.bioinf.ms.middleware.model.annotations;

import de.unijena.bioinf.ChemistryBase.ms.AnnotatedPeak;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.ft.*;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.Getter;
import lombok.Setter;

import java.util.*;
import java.util.stream.Collectors;
//todo create tree serialization that does not store every node multiple times
@Getter
@Setter
public class FragmentationTree {
    protected List<FragmentNode> fragments;
    protected List<LossEdge> losses;
    Double treeScore;
    FragmentNode root;

    public static FragmentationTree fromFtree(FTree sourceTree) {
        final FragmentationTree tree = new FragmentationTree();
        Map<Integer, FragmentNode> fragmentNodes = new HashMap<>();
        List<LossEdge> lossEdges = new ArrayList<>();

        tree.setTreeScore(sourceTree.getTreeWeight());

        for (Fragment f : sourceTree.getFragments()) {
            final FragmentNode fn = new FragmentNode();
            fn.setId(f.getVertexId());
            fn.setMolecularFormula(f.getFormula().toString());
            fn.setIonType(f.getIonization().toString());

            {
                final FragmentAnnotation<Peak> peakInfo = sourceTree.getFragmentAnnotationOrThrow(Peak.class);
                fn.setMz(peakInfo.get(f).getMass());
                fn.setIntensity(peakInfo.get(f).getIntensity());
            }

            {
                final FragmentAnnotation<AnnotatedPeak> anoPeak = sourceTree.getFragmentAnnotationOrThrow(AnnotatedPeak.class);
                if (anoPeak.get(f).isMeasured()) {
                    Deviation dev = sourceTree.getMassError(f);
                    fn.setMassDeviationDa(dev.getAbsolute());
                    fn.setMassDeviationPpm(dev.getPpm());
                }
            }

            {
                final FragmentAnnotation<Score> scores = sourceTree.getFragmentAnnotationOrThrow(Score.class);
                fn.setScore(scores.get(f).sum());
            }

            fragmentNodes.put(fn.getId(), fn);
        }

        tree.setFragments(fragmentNodes.values().stream().sorted(Comparator.comparing(FragmentNode::getId)).collect(Collectors.toList()));
        tree.setRoot(fragmentNodes.get(0));

        for (Loss l : sourceTree.losses()) {
            LossEdge loss = new LossEdge();
            loss.setSourceFragment(fragmentNodes.get(l.getSource().getVertexId()));
            loss.setTargetFragment(fragmentNodes.get(l.getTarget().getVertexId()));
            loss.setMolecularFormula(l.getFormula().toString());
            loss.setScore(l.getWeight());
            lossEdges.add(loss);
        }
        tree.setLosses(lossEdges);

        return tree;
    }
}