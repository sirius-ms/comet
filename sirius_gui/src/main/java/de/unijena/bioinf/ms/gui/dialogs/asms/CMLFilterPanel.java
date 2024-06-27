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
import java.util.Optional;

/**
 * Compound filter options for measured combinatorial molecule libraries (especially for affinity selection experiments).
 */
public class CMLFilterPanel extends JPanel {

    private final CompoundFilterModel filterModel;
    private final FileChooserPanel bbFileSelectionPanel, outputFileSelectionPanel;
    private final JTextField scaffoldMolFormulaField;
    private final JSpinner minPeaksSpinner, numTopPeaksSpinner, ms1DevSpinner, ms2DevSpinner, numHydrogenShiftsSpinner;
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
                this.bbFileSelectionPanel = new FileChooserPanel(Optional.ofNullable(cmlFilterModel.getPathToBBFile()).orElse(""));
                this.scaffoldMolFormulaField = new JTextField(Optional.ofNullable(cmlFilterModel.getScaffoldMf()).orElse(""));

                Box libSelBox = Box.createHorizontalBox();
                libSelBox.add(new TwoColumnPanel("Building blocks:", bbFileSelectionPanel));
                libSelBox.add(Box.createHorizontalGlue());
                libSelBox.add(new TwoColumnPanel("Scaffold formula:", scaffoldMolFormulaField));
                params.add(libSelBox);
            }

            // MS1 deviation:
            {
                this.ms1DevSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getMs1Deviation(), 0, Double.POSITIVE_INFINITY, 1d));
                params.addNamed("MS1 mass accuracy (ppm):", this.ms1DevSpinner);
            }
        }

        // Peak matching filter params:
        {
            // Params for enabling this filter
            {
                this.peakMatchingFilterCheckBox = new JCheckBox("Enable peak matching filter", cmlFilterModel.isPeakMatchingFilterEnabled());
                this.peakMatchingFilterCheckBox.addChangeListener(new PeakMatchingFilterCheckBoxListener());

                params.add(Box.createVerticalStrut(5));
                Box box = Box.createHorizontalBox();
                box.add(this.peakMatchingFilterCheckBox);
                box.add(Box.createHorizontalStrut(25));
                params.add(box);
            }

            // Params for the minimal and maximal number of peaks to consider:
            {
                this.minPeaksSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getMinMatchingPeaks(), 0, cmlFilterModel.getNumTopPeaks(), 1));
                this.numTopPeaksSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getNumTopPeaks(), 0, Integer.MAX_VALUE, 1));
                this.numTopPeaksSpinner.addChangeListener(new NumTopPeaksListener());

                Box numMatchingPeaksBox = Box.createHorizontalBox();
                numMatchingPeaksBox.add(new TwoColumnPanel("Minimum number of matching peaks:", this.minPeaksSpinner));
                numMatchingPeaksBox.add(Box.createHorizontalGlue());
                numMatchingPeaksBox.add(new TwoColumnPanel("Number of considered peaks:", this.numTopPeaksSpinner));
                params.add(numMatchingPeaksBox);
            }

            // Parameters for the mass accuracy and the location of the output file:
            {
                this.numHydrogenShiftsSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getNumAllowedHydrogenShifts(), 0, Integer.MAX_VALUE, 1));
                this.ms2DevSpinner = new JSpinner(new SpinnerNumberModel(cmlFilterModel.getMs2Deviation(), 0, Double.POSITIVE_INFINITY, 1d));
                this.outputFileSelectionPanel = new FileChooserPanel();

                Box box = Box.createHorizontalBox();
                box.add(new TwoColumnPanel("Number of allowed hydrogen shifts:", this.numHydrogenShiftsSpinner));
                box.add(Box.createHorizontalGlue());
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
        this.numHydrogenShiftsSpinner.setEnabled(isEnabled);
        this.outputFileSelectionPanel.setEnabled(isEnabled);
    }

    public void applyToModel(@NotNull CompoundFilterModel filterModel){
        final String pathToBBFile = Optional.ofNullable(this.bbFileSelectionPanel.getFilePath()).map(String::strip).orElse(null);
        final String outputPath = Optional.ofNullable(this.outputFileSelectionPanel.getFilePath()).map(String::strip).orElse(null);
        final String scaffoldMf = this.getStringValue(this.scaffoldMolFormulaField); // return null if text is null or consists of only whitespaces
        final boolean isPeakMatchingFilterEnabled = this.peakMatchingFilterCheckBox.isSelected();
        final int minMatchingPeaks = this.getIntValue(this.minPeaksSpinner);
        final int numTopPeaks = this.getIntValue(this.numTopPeaksSpinner);
        final double ms1Dev = this.getDoubleValue(this.ms1DevSpinner);
        final double ms2Dev = this.getDoubleValue(this.ms2DevSpinner);
        final int numHydrogenShifts = this.getIntValue(this.numHydrogenShiftsSpinner);

        CMLFilterModelOptions cmlFilterOptions = new CMLFilterModelOptions(pathToBBFile, scaffoldMf, outputPath, minMatchingPeaks,
                numTopPeaks, numHydrogenShifts, ms1Dev, ms2Dev, isPeakMatchingFilterEnabled);
        filterModel.setCMLFilterOptions(cmlFilterOptions);
    }

    public void reset(){
        final CMLFilterModelOptions defaultOptions = CMLFilterModelOptions.disabled();
        this.bbFileSelectionPanel.field.setText(defaultOptions.getPathToBBFile());
        this.outputFileSelectionPanel.field.setText(defaultOptions.getMatchedPeaksOutputFilePath());
        this.scaffoldMolFormulaField.setText(defaultOptions.getScaffoldMf());
        this.minPeaksSpinner.setValue(defaultOptions.getMinMatchingPeaks());
        this.numTopPeaksSpinner.setValue(defaultOptions.getNumTopPeaks());
        this.numHydrogenShiftsSpinner.setValue(defaultOptions.getNumAllowedHydrogenShifts());
        this.ms1DevSpinner.setValue(defaultOptions.getMs1Deviation());
        this.ms2DevSpinner.setValue(defaultOptions.getMs2Deviation());
        this.peakMatchingFilterCheckBox.setSelected(defaultOptions.isPeakMatchingFilterEnabled());
    }

    private double getDoubleValue(JSpinner spinner){
        return ((SpinnerNumberModel) spinner.getModel()).getNumber().doubleValue();
    }

    private int getIntValue(JSpinner spinner){
        return ((SpinnerNumberModel) spinner.getModel()).getNumber().intValue();
    }

    private String getStringValue(JTextField textField){
        final String txt = textField.getText();
        if(txt == null || txt.isBlank())
            return null;
        return txt.strip();
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
