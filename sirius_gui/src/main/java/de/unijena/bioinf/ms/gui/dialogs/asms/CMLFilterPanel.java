package de.unijena.bioinf.ms.gui.dialogs.asms;

import de.unijena.bioinf.ms.frontend.io.FileChooserPanel;
import de.unijena.bioinf.ms.gui.utils.CompoundFilterModel;
import de.unijena.bioinf.ms.gui.utils.TwoColumnPanel;
import de.unijena.bioinf.ms.gui.utils.asms.CMLFilterModelOptions;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;

/**
 * Compound filter options for measured combinatorial molecule libraries (especially for affinity selection experiments).
 */
public class CMLFilterPanel extends JPanel {

    private CompoundFilterModel filterModel;
    private FileChooserPanel bbFileSelectionPanel, matchedPeaksOutputFile;
    private JTextField scaffoldMolFormulaField;
    private JSpinner minPeaksSpinner, numTopPeaksSpinner, ms1DevSpinner, ms2DevSpinner;

    public CMLFilterPanel(CompoundFilterModel filterModel){
        this.filterModel = filterModel;
        CMLFilterModelOptions cmlFilterModel = this.filterModel.getCmlFilterOptions();

        TwoColumnPanel params = new TwoColumnPanel();
        this.add(params, BorderLayout.CENTER);

        // Library params: BB-File + Scaffold molecular formula
        this.bbFileSelectionPanel = new FileChooserPanel();
        this.matchedPeaksOutputFile = new FileChooserPanel();
        this.scaffoldMolFormulaField = new JTextField();

        Box libSelBox = Box.createHorizontalBox();
        libSelBox.add(new TwoColumnPanel("Building blocks:", bbFileSelectionPanel));
        libSelBox.add(Box.createHorizontalGlue());
        libSelBox.add(new TwoColumnPanel("Scaffold formula:", scaffoldMolFormulaField));
        params.add(libSelBox);
        params.addNamed("Output location:", matchedPeaksOutputFile);

        // Peak matching filter params:
        this.minPeaksSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentMinMatchingPeaks(), 0, cmlFilterModel.getCurrentNumTopPeaks(), 1));
        this.numTopPeaksSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentNumTopPeaks(), 0, Integer.MAX_VALUE, 1));
        this.numTopPeaksSpinner.addChangeListener(new NumTopPeaksListener());

        Box numMatchingPeaksBox = Box.createHorizontalBox();
        numMatchingPeaksBox.add(new TwoColumnPanel("Minimum number of matching peaks:", this.minPeaksSpinner));
        numMatchingPeaksBox.add(new TwoColumnPanel("Number of considered peaks:", this.numTopPeaksSpinner));
        params.add(numMatchingPeaksBox);

        this.ms1DevSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentMs1Deviation(), 0,Double.POSITIVE_INFINITY, 1d));
        this.ms2DevSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentMs2Deviation(), 0, Double.POSITIVE_INFINITY, 1d));

        params.addNamed("MS1 mass accuracy (ppm)", this.ms1DevSpinner);
        params.addNamed("MS2 mass accuracy (ppm)", this.ms2DevSpinner);

        this.add(params);
    }


    private class NumTopPeaksListener implements ChangeListener {

        @Override
        public void stateChanged(ChangeEvent e) {
            int currentMinPeaksValue = (int) minPeaksSpinner.getValue();
            int currentNumTopPeaksValue = (int) numTopPeaksSpinner.getValue();
            if(currentMinPeaksValue > currentNumTopPeaksValue){
                minPeaksSpinner.setModel(new SpinnerNumberModel(currentNumTopPeaksValue, 0, currentNumTopPeaksValue, 1));
            }
        }
    }


}
