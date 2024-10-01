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

import ca.odell.glazedlists.swing.DefaultEventListModel;
import de.unijena.bioinf.ms.gui.SiriusGui;
import de.unijena.bioinf.ms.gui.actions.SiriusActions;
import de.unijena.bioinf.ms.gui.utils.JListDropImage;
import de.unijena.bioinf.projectspace.InstanceBean;

import javax.swing.*;
import java.awt.event.*;

/**
 * @author Markus Fleischauer (markus.fleischauer@gmail.com)
 */
public class CompoundListView extends JScrollPane {

    final CompoundList sourceList;
    final JListDropImage<InstanceBean> compoundListView;
    final JPopupMenu expPopMenu;

    public CompoundListView(SiriusGui gui, CompoundList sourceList) {
        super(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        this.sourceList = sourceList;
        //todo move texfield and filter funktion here
        compoundListView = new JListDropImage<>(new DefaultEventListModel<>(sourceList.compoundList));
        compoundListView.setSelectionModel(sourceList.compountListSelectionModel);
        compoundListView.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        compoundListView.setCellRenderer(new CompoundCellRenderer(gui));
        expPopMenu = new CompoundContextMenu(gui);
        compoundListView.addMouseListener(new MouseListener() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    // Double-click detected
                    int index = compoundListView.locationToIndex(e.getPoint());
                    compoundListView.setSelectedIndex(index);
                    SiriusActions.COMPUTE.getInstance(gui, true)
                            .actionPerformed(new ActionEvent(compoundListView, 123, SiriusActions.COMPUTE.name()));
                }
            }

            @Override
            public void mousePressed(MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e)) {
                    int indx = compoundListView.locationToIndex(e.getPoint());
                    boolean select = true;
                    for (int i : compoundListView.getSelectedIndices()) {
                        if (indx == i) {
                            select = false;
                            break;
                        }
                    }
                    if (select) {
                        compoundListView.setSelectedIndex(indx);
                    }

                    if (e.isPopupTrigger()) {
                        expPopMenu.show(e.getComponent(), e.getX(), e.getY());
                    }
                }
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                if (e.isPopupTrigger()) {
                    expPopMenu.show(e.getComponent(), e.getX(), e.getY());
                }
            }

            @Override
            public void mouseEntered(MouseEvent e) {
            }

            @Override
            public void mouseExited(MouseEvent e) {
            }
        });


        setViewportView(compoundListView);

        KeyStroke enterKey = KeyStroke.getKeyStroke("ENTER");
        compoundListView.getInputMap().put(enterKey, SiriusActions.COMPUTE.name());


        KeyStroke delKey = KeyStroke.getKeyStroke("DELETE");
        compoundListView.getInputMap().put(delKey, SiriusActions.DELETE_EXP.name());


        ActionMap actionMap = compoundListView.getActionMap();

        // Define and register the compute action (for Enter key)
        actionMap.put(SiriusActions.COMPUTE.name(), new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                SiriusActions.COMPUTE.getInstance(gui, true)
                        .actionPerformed(new ActionEvent(compoundListView, ActionEvent.ACTION_PERFORMED, SiriusActions.COMPUTE.name()));
            }
        });

        // Define and register the delete action (for Delete key)
        actionMap.put(SiriusActions.DELETE_EXP.name(), new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                SiriusActions.DELETE_EXP.getInstance(gui, true)
                        .actionPerformed(new ActionEvent(compoundListView, ActionEvent.ACTION_PERFORMED, SiriusActions.DELETE_EXP.name()));
            }
        });
    }

    public void ensureIndexIsVisible(int index) {
        compoundListView.ensureIndexIsVisible(index);
    }

    public JPopupMenu getExpPopMenu() {
        return expPopMenu;
    }

}
