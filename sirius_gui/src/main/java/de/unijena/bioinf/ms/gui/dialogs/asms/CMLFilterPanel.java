package de.unijena.bioinf.ms.gui.dialogs.asms;

import de.unijena.bioinf.ms.frontend.io.FileChooserPanel;
import de.unijena.bioinf.ms.gui.utils.CompoundFilterModel;
import de.unijena.bioinf.ms.gui.utils.TwoColumnPanel;
import de.unijena.bioinf.ms.gui.utils.asms.CMLFilterModelOptions;
import org.jetbrains.annotations.NotNull;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;

/**
 * Compound filter options for measured combinatorial molecule libraries (especially for affinity selection experiments).
 */
public class CMLFilterPanel extends JPanel {

    private final CompoundFilterModel filterModel;
    private final FileChooserPanel bbFileSelectionPanel, outputFileSelectionPanel;
    private final JTextField scaffoldMolFormulaField;
    private final JSpinner minPeaksSpinner, numTopPeaksSpinner, ms1DevSpinner, ms2DevSpinner;
    private final JCheckBox peakMatchingFilterCheckBox;

    public CMLFilterPanel(CompoundFilterModel filterModel){
        this.filterModel = filterModel;
        CMLFilterModelOptions cmlFilterModel = this.filterModel.getCmlFilterOptions();

        TwoColumnPanel params = new TwoColumnPanel();
        this.add(params, BorderLayout.CENTER);

        // Ms1 candidate filter params:
        {
            // Library params: BB file + scaffold molecular formula
            {
                this.bbFileSelectionPanel = new FileChooserPanel();
                this.scaffoldMolFormulaField = new JTextField();

                Box libSelBox = Box.createHorizontalBox();
                libSelBox.add(new TwoColumnPanel("Building blocks:", bbFileSelectionPanel));
                libSelBox.add(Box.createHorizontalGlue());
                libSelBox.add(new TwoColumnPanel("Scaffold formula:", scaffoldMolFormulaField));
                params.add(libSelBox);
            }

            // MS1 deviation:
            {
                this.ms1DevSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentMs1Deviation(), 0, Double.POSITIVE_INFINITY, 1d));
                params.addNamed("MS1 mass accuracy (ppm):", this.ms1DevSpinner);
            }
        }

        // Peak matching filter params:
        {
            // Params for enabling this filter
            {
                this.peakMatchingFilterCheckBox = new JCheckBox("Enable peak matching filter");
                this.peakMatchingFilterCheckBox.addChangeListener(new PeakMatchingFilterCheckBoxListener());

                params.add(Box.createVerticalStrut(5));
                Box box = Box.createHorizontalBox();
                box.add(this.peakMatchingFilterCheckBox);
                box.add(Box.createHorizontalStrut(25));
                params.add(box);
            }

            // Params for the minimal and maximal number of peaks to consider:
            {
                this.minPeaksSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentMinMatchingPeaks(), 0, cmlFilterModel.getCurrentNumTopPeaks(), 1));
                this.numTopPeaksSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentNumTopPeaks(), 0, Integer.MAX_VALUE, 1));
                this.numTopPeaksSpinner.addChangeListener(new NumTopPeaksListener());

                Box numMatchingPeaksBox = Box.createHorizontalBox();
                numMatchingPeaksBox.add(new TwoColumnPanel("Minimum number of matching peaks:", this.minPeaksSpinner));
                numMatchingPeaksBox.add(Box.createHorizontalGlue());
                numMatchingPeaksBox.add(new TwoColumnPanel("Number of considered peaks:", this.numTopPeaksSpinner));
                params.add(numMatchingPeaksBox);
            }

            // Parameters for the mass accuracy and the location of the output file:
            {
                this.ms2DevSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getCurrentMs2Deviation(), 0, Double.POSITIVE_INFINITY, 1d));
                this.outputFileSelectionPanel = new FileChooserPanel();

                Box box = Box.createHorizontalBox();
                box.add(new TwoColumnPanel("MS2 mass accuracy (ppm):", this.ms2DevSpinner));
                box.add(Box.createHorizontalGlue());
                box.add(new TwoColumnPanel("Output location:", this.outputFileSelectionPanel));
                params.add(box);
            }
        }
        this.setEnablePeakMatchingFilter(false);
        this.add(params);
    }

    private void setEnablePeakMatchingFilter(boolean isEnabled){
        this.minPeaksSpinner.setEnabled(isEnabled);
        this.numTopPeaksSpinner.setEnabled(isEnabled);
        this.ms2DevSpinner.setEnabled(isEnabled);
        this.outputFileSelectionPanel.setEnabled(isEnabled);
    }

    /* todo
    public void applyToModel(@NotNull CompoundFilterModel filterModel){
        String pathToBBFile = this.bbFileSelectionPanel.getFilePath();
        String outputPath = this.outputFileSelectionPanel.getFilePath();
        String scaffoldMf = this.scaffoldMolFormulaField.getText();
        int numTopPeaks = this.getIntValue(this.numTopPeaksSpinner);
        double ms1Dev = this.getDoubleValue(this.ms1DevSpinner);
        double ms2Dev = this.getDoubleValue(this.ms2DevSpinner);

        CMLFilterModelOptions cmlFilterOptions = new CMLFilterModelOptions(pathToBBFile, scaffoldMf, outputPath,
                this.getIntValue(this.minPeaksSpinner), this.getIntValue(this.numTopPeaksSpinner), )


    }

     */


    private double getDoubleValue(JSpinner spinner){
        return ((SpinnerNumberModel) spinner.getModel()).getNumber().doubleValue();
    }

    private int getIntValue(JSpinner spinner){
        return ((SpinnerNumberModel) spinner.getModel()).getNumber().intValue();
    }

    private class NumTopPeaksListener implements ChangeListener {

        @Override
        public void stateChanged(ChangeEvent e) {
            int currentMinPeaksValue = ((SpinnerNumberModel) minPeaksSpinner.getModel()).getNumber().intValue();
            int currentNumTopPeaksValue = ((SpinnerNumberModel) numTopPeaksSpinner.getModel()).getNumber().intValue();
            if(currentMinPeaksValue > currentNumTopPeaksValue) minPeaksSpinner.getModel().setValue(currentNumTopPeaksValue);
            ((SpinnerNumberModel) minPeaksSpinner.getModel()).setMaximum(currentNumTopPeaksValue);
        }
    }

    private class PeakMatchingFilterCheckBoxListener implements ChangeListener {

        @Override
        public void stateChanged(ChangeEvent e) {
            final boolean isFilterEnabled = peakMatchingFilterCheckBox.isSelected();
            setEnablePeakMatchingFilter(isFilterEnabled);
        }
    }


}
