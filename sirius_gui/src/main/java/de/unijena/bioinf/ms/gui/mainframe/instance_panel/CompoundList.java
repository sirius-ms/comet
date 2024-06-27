/*
 *  This file is part of the SIRIUS Software for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer, Marvin Meusel and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schiller University.
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Affero General Public License
 *  as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License along with SIRIUS.  If not, see <https://www.gnu.org/licenses/agpl-3.0.txt>
 */

package de.unijena.bioinf.ms.gui.mainframe.instance_panel;

import ca.odell.glazedlists.*;
import ca.odell.glazedlists.event.ListEvent;
import ca.odell.glazedlists.matchers.CompositeMatcherEditor;
import ca.odell.glazedlists.matchers.MatcherEditor;
import ca.odell.glazedlists.swing.DefaultEventSelectionModel;
import ca.odell.glazedlists.swing.GlazedListsSwing;
import ca.odell.glazedlists.swing.TextComponentMatcherEditor;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.datastructures.BuildingBlock;
import de.unijena.bioinf.datastructures.CMLCandidates;
import de.unijena.bioinf.datastructures.OrderedCombinatorialMoleculeLibrary;
import de.unijena.bioinf.datastructures.Scaffold;
import de.unijena.bioinf.io.BuildingBlockReader;
import de.unijena.bioinf.ms.gui.SiriusGui;
import de.unijena.bioinf.ms.gui.dialogs.CompoundFilterOptionsDialog;
import de.unijena.bioinf.ms.gui.utils.*;
import de.unijena.bioinf.ms.gui.utils.matchers.BackgroundJJobMatcheEditor;
import de.unijena.bioinf.projectspace.GuiProjectManager;
import de.unijena.bioinf.projectspace.InstanceBean;
import lombok.Getter;
import org.jetbrains.annotations.NotNull;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * This is the main List of the SIRIUS UI.
 * It shows the main Instances (former Compounds or Experiments)
 * It is usually a singleton and backed by the INSTANCE_LIST of the  {@link GuiProjectManager}
 *
 * @author Markus Fleischauer (markus.fleischauer@gmail.com)
 */
public class CompoundList {

    final PlaceholderTextField searchField;
    final JButton openFilterPanelButton;
    final CompoundFilterModel compoundFilterModel;
    final ObservableElementList<InstanceBean> obsevableScource;
    final SortedList<InstanceBean> sortedSource;
    final EventList<InstanceBean> compoundList;
    final DefaultEventSelectionModel<InstanceBean> compountListSelectionModel;
    final BackgroundJJobMatcheEditor<InstanceBean> backgroundFilterMatcher;
    final private MatcherEditorWithOptionalInvert<InstanceBean> compoundListMatchEditor;

    @Getter
    private OrderedCombinatorialMoleculeLibrary cmlLibrary;

    private final Queue<ExperimentListChangeListener> listeners = new ConcurrentLinkedQueue<>();

    private final Color defaultOpenFilterPanelButtonColor;

    @Getter
    private @NotNull SiriusGui gui;
    public CompoundList(@NotNull SiriusGui gui) {
        this.gui = gui;
        searchField = new PlaceholderTextField();
        searchField.setPlaceholder("Type and hit enter to search");
        searchField.setToolTipText("Type text to perform a full text search on the data below. Hit enter to start searching.");

        obsevableScource = new ObservableElementList<>(gui.getProjectManager().INSTANCE_LIST, GlazedLists.beanConnector(InstanceBean.class));
        sortedSource = new SortedList<>(obsevableScource, Comparator.comparing(InstanceBean::getRTOrMissing));

        //filters
        BasicEventList<MatcherEditor<InstanceBean>> listOfFilters = new BasicEventList<>();
        //text filter
        listOfFilters.add(new TextComponentMatcherEditor<>(searchField, (baseList, element) -> {
            baseList.add(element.getGUIName());
            baseList.add(element.getIonType().toString());
            baseList.add(String.valueOf(element.getIonMass()));
        }, false));

        //additional filter based on specific parameters
        compoundFilterModel = new CompoundFilterModel(this);
        listOfFilters.add(new CompoundFilterMatcherEditor(new CompoundFilterMatcher(gui.getProperties(), compoundFilterModel)));
        //combined filters
        CompositeMatcherEditor<InstanceBean> compositeMatcherEditor = new CompositeMatcherEditor<>(listOfFilters);
        compositeMatcherEditor.setMode(CompositeMatcherEditor.AND);

        compoundListMatchEditor = new MatcherEditorWithOptionalInvert<>(compositeMatcherEditor);
        backgroundFilterMatcher = new BackgroundJJobMatcheEditor<>(compoundListMatchEditor);
        FilterList<InstanceBean> filterList = new FilterList<>(sortedSource, backgroundFilterMatcher);
        compoundList = GlazedListsSwing.swingThreadProxyList(filterList);

        //filter dialog
        openFilterPanelButton = new JButton("...");
        defaultOpenFilterPanelButtonColor = openFilterPanelButton.getBackground();

        openFilterPanelButton.addActionListener(e -> new CompoundFilterOptionsDialog(gui, searchField, compoundFilterModel, this));
        compositeMatcherEditor.addMatcherEditorListener(evt -> colorByActiveFilter(openFilterPanelButton, compoundFilterModel));

        compountListSelectionModel = new DefaultEventSelectionModel<>(compoundList);

        compountListSelectionModel.addListSelectionListener(e -> {
            if (!e.getValueIsAdjusting()) {
                compountListSelectionModel.getDeselected().forEach(InstanceBean::disableProjectSpaceListener);
                compountListSelectionModel.getSelected().forEach(InstanceBean::enableProjectSpaceListener);
                notifyListenerSelectionChange();
            }
        });
        compoundList.addListEventListener(this::notifyListenerDataChange);
    }

    private void colorByActiveFilter(JButton openFilterPanelButton, CompoundFilterModel compoundFilterModel) {
        //is any filtering option active (despite the text filter which is visible all the time)
        if (isFilterInverted()){
            openFilterPanelButton.setBackground(new Color(235, 94, 85));
        } else if (compoundFilterModel.isActive() || !searchField.getText().isEmpty()){
            openFilterPanelButton.setBackground(new Color(49, 153, 187));
        }else {
            openFilterPanelButton.setBackground(defaultOpenFilterPanelButtonColor);
        }
    }

    public void orderBy(@NotNull final Comparator<InstanceBean> comp) {
        sortedSource.setComparator(comp);
    }

    public boolean isFilterInverted() {
        return compoundListMatchEditor.isInverted();
    }

    public void toggleInvertFilter() {
        compoundListMatchEditor.setInverted(!compoundListMatchEditor.isInverted());
    }

    public void resetFilter() {
        //filtering consists of the text filter, the filter model and the possible inversion using the MatcherEditor
        compoundFilterModel.resetFilter();
        compoundListMatchEditor.setInverted(false);
        searchField.setText("");
        searchField.postActionEvent();
        colorByActiveFilter(openFilterPanelButton, compoundFilterModel);
    }

    private void notifyListenerDataChange(ListEvent<InstanceBean> event) {
        for (ExperimentListChangeListener l : listeners) {
            event.reset();//this is hell important to reset the iterator
            l.listChanged(event, compountListSelectionModel, sortedSource.size());
        }
    }

    private void notifyListenerSelectionChange() {
        for (ExperimentListChangeListener l : listeners) {
            l.listSelectionChanged(compountListSelectionModel, sortedSource.size());
        }
    }

    //API methods
    public void addChangeListener(ExperimentListChangeListener l) {
        listeners.add(l);
    }

    public void removeChangeListener(ExperimentListChangeListener l) {
        listeners.remove(l);
    }

    public DefaultEventSelectionModel<InstanceBean> getCompoundListSelectionModel() {
        return compountListSelectionModel;
    }

    public EventList<InstanceBean> getCompoundList() {
        return compoundList;
    }

    public void initCmlLibraryAndRemoveCMLAnnotations(String pathToBBFile, String scaffoldMf){
        // 1. Initialize the combinatorial molecule library:
        if (pathToBBFile != null) {
            try {
                final File bbFile = new File(pathToBBFile);
                final BuildingBlock[][] buildingBlocks = BuildingBlockReader.readBuildingBlocks(bbFile);
                final Scaffold scaffold = new Scaffold(scaffoldMf != null ? MolecularFormula.parse(scaffoldMf) : MolecularFormula.emptyFormula(), null);
                cmlLibrary = new OrderedCombinatorialMoleculeLibrary(bbFile.getName().replaceFirst(".csv", ""), buildingBlocks, scaffold);
            } catch (UnknownElementException | IOException e) {
                e.printStackTrace();
                cmlLibrary = OrderedCombinatorialMoleculeLibrary.emptyLibrary();
            }
        } else {
            cmlLibrary = OrderedCombinatorialMoleculeLibrary.emptyLibrary();
        }
         // 2. Delete previous Ms2Experiment annotation of CMLCandidates.class:
        this.removeCMLAnnotations();
    }

    public void removeCMLAnnotations(){
        compoundList.forEach(item -> {
            final Ms2Experiment exp = item.asMs2Experiment();
            exp.removeAnnotation(CMLCandidates.class);
        });
    }
}
