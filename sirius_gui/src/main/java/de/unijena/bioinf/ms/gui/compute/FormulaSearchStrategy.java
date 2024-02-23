package de.unijena.bioinf.ms.gui.compute;

import de.unijena.bioinf.ChemistryBase.chem.ChemicalAlphabet;
import de.unijena.bioinf.ChemistryBase.chem.Element;
import de.unijena.bioinf.ChemistryBase.chem.FormulaConstraints;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.ChemistryBase.ms.ft.model.FormulaSettings;
import de.unijena.bioinf.ChemistryBase.utils.DescriptiveOptions;
import de.unijena.bioinf.chemdb.annotations.FormulaSearchDB;
import de.unijena.bioinf.chemdb.custom.CustomDataSources;
import de.unijena.bioinf.ms.frontend.core.ApplicationCore;
import de.unijena.bioinf.ms.frontend.subtools.sirius.SiriusOptions;
import de.unijena.bioinf.ms.gui.dialogs.ElementSelectionDialog;
import de.unijena.bioinf.ms.gui.utils.GuiUtils;
import de.unijena.bioinf.ms.gui.utils.TextHeaderBoxPanel;
import de.unijena.bioinf.ms.gui.utils.TwoColumnPanel;
import de.unijena.bioinf.ms.gui.utils.jCheckboxList.JCheckboxListPanel;
import de.unijena.bioinf.ms.properties.PropertyManager;
import de.unijena.bioinf.projectspace.InstanceBean;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;

public class FormulaSearchStrategy extends ConfigPanel {
    public enum Strategy implements DescriptiveOptions {
        DEFAULT("Default strategy"),
        DE_NOVO("Denovo strategy"),
        DATABASE("Database stragegy");

        private final String description;

        Strategy(String description) {
            this.description = description;
        }

        @Override
        public String getDescription() {
            return description;
        }
    }

    public enum ElementAlphabetStrategy implements DescriptiveOptions {
        DE_NOVO_ONLY("Use set of elements for de novo generation only."),
        BOTH("Use set of elements for de novo generation and filter of bottom up search.");

        private final String description;

        ElementAlphabetStrategy(String description) {
            this.description = description;
        }

        @Override
        public String getDescription() {
            return description;
        }
    }

    protected Strategy strategy;

    protected final Dialog owner;
    protected final List<InstanceBean> ecs;
    protected final boolean isMs2;
    protected final boolean isBatchDialog;

    protected  JCheckboxListPanel<CustomDataSources.Source> searchDBList;

    private JPanel elementFilterPanel;

    /**
     * Map of strategy-specific UI components for showing/hiding when changing the strategy
     */
    private final Map<Strategy, List<Component>> strategyComponents;

    public FormulaSearchStrategy(Dialog owner, List<InstanceBean> ecs, boolean isMs2, boolean isBatchDialog, ParameterBinding parameterBindings) {
        super(parameterBindings);
        this.owner = owner;
        this.ecs = ecs;
        this.isMs2 = isMs2;
        this.isBatchDialog = isBatchDialog;

        strategyComponents = new HashMap<>();
        strategyComponents.put(Strategy.DEFAULT, new ArrayList<>());
        strategyComponents.put(Strategy.DE_NOVO, new ArrayList<>());
        strategyComponents.put(Strategy.DATABASE, new ArrayList<>());

        createPanel();
    }

    public JCheckboxListPanel<CustomDataSources.Source> getSearchDBList() {
        return searchDBList;
    }

    private void createPanel() {
        setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));

        final JPanel formulaSearchStrategySelection = new JPanel();
        formulaSearchStrategySelection.setLayout(new BoxLayout(formulaSearchStrategySelection, BoxLayout.PAGE_AXIS));
        formulaSearchStrategySelection.setBorder(BorderFactory.createEmptyBorder(0, GuiUtils.LARGE_GAP, 0, 0));
        JComboBox<FormulaSearchStrategy.Strategy> strategyBox =  GuiUtils.makeParameterComboBoxFromDescriptiveValues(FormulaSearchStrategy.Strategy.values());
        formulaSearchStrategySelection.add(new TextHeaderBoxPanel("Molecular formula generation", strategyBox));

        add(formulaSearchStrategySelection);
        add(Box.createRigidArea(new Dimension(0, GuiUtils.MEDIUM_GAP)));

        JPanel strategyCardContainer = new JPanel();
        strategyCardContainer.setBorder(BorderFactory.createEmptyBorder(0, GuiUtils.LARGE_GAP, 0, 0));
        strategyCardContainer.setLayout(new BoxLayout(strategyCardContainer, BoxLayout.PAGE_AXIS));

        strategy = (Strategy) strategyBox.getSelectedItem();
        elementFilterPanel = createElementFilterPanel();

        JPanel defaultStrategyParameters = createDefaultStrategyParameters();
        JPanel databaseStrategyParameters = createDatabaseStrategyParameters();

        strategyComponents.get(Strategy.DEFAULT).add(defaultStrategyParameters);
        strategyComponents.get(Strategy.DATABASE).add(databaseStrategyParameters);

        strategyCardContainer.add(defaultStrategyParameters);
        strategyCardContainer.add(databaseStrategyParameters);

        strategyCardContainer.add(createElementFilterPanel());

        add(strategyCardContainer);

        showStrategy(strategy);

        strategyBox.addItemListener(e -> {
            if (e.getStateChange() != ItemEvent.SELECTED) {
                return;
            }
            final Strategy source = (Strategy) e.getItem();
            showStrategy(source);
        });
    }

    private void showStrategy(Strategy strategy) {
        strategyComponents.forEach((s, lst) -> lst.forEach(c -> c.setVisible(s.equals(strategy))));
    }

    private JPanel createDefaultStrategyParameters() {
        JPanel card = new JPanel();
        card.setLayout(new BoxLayout(card, BoxLayout.PAGE_AXIS));
        JPanel parameterPanel = applyDefaultLayout(new JPanel());

        final TwoColumnPanel busOptions = new TwoColumnPanel();
        JComboBox<SiriusOptions.BottomUpSearchOptions> bottomUpSearchSelector = new JComboBox<>();
        List<SiriusOptions.BottomUpSearchOptions> settings = new ArrayList<>(EnumSet.allOf(SiriusOptions.BottomUpSearchOptions.class));
        if (strategy == Strategy.DE_NOVO) settings.remove(SiriusOptions.BottomUpSearchOptions.BOTTOM_UP_ONLY);
        if (strategy == Strategy.DEFAULT) settings.remove(SiriusOptions.BottomUpSearchOptions.DISABLED);//this is not a contradiction by default, but we have the separate Strategy.DE_NOVO for that
        settings.forEach(bottomUpSearchSelector::addItem);
        busOptions.addNamed("Bottom up search", bottomUpSearchSelector);

        JSpinner bottomUpSearchOnly;
        JCheckBox bottomUpSearchEnabled = new JCheckBox();
        bottomUpSearchOnly = makeIntParameterSpinner("FormulaSearchSettings.disableDeNovoAboveMass", 0, Integer.MAX_VALUE, 5);
        bottomUpSearchEnabled.setEnabled(false);
        bottomUpSearchOnly.setEnabled(false);

        bottomUpSearchEnabled.addActionListener(e -> {
            //enable of disable buttom up search. Not setting specific mz threshold
            setBottomUpSearchMass(bottomUpSearchEnabled);
        });
        bottomUpSearchEnabled.setSelected(true);

        bottomUpSearchSelector.addItemListener(e -> {
            if (e.getStateChange() != ItemEvent.SELECTED) {
                return;
            }
            final SiriusOptions.BottomUpSearchOptions source = (SiriusOptions.BottomUpSearchOptions) e.getItem();
            switch (source) {
                case BOTTOM_UP_ONLY -> {
                    bottomUpSearchEnabled.setEnabled(false);
                    bottomUpSearchOnly.setEnabled(false);
                    bottomUpSearchEnabled.setSelected(true);
                    bottomUpSearchOnly.setValue(0);
                }
                case DISABLED -> {
                    bottomUpSearchEnabled.setEnabled(false);
                    bottomUpSearchOnly.setEnabled(false);
                    bottomUpSearchEnabled.setSelected(false);
                    bottomUpSearchOnly.setValue(Double.POSITIVE_INFINITY);
                }
                case CUSTOM -> {
                    bottomUpSearchEnabled.setEnabled(true);
                    bottomUpSearchOnly.setEnabled(true);
                    bottomUpSearchEnabled.setSelected(true);
                    bottomUpSearchOnly.setValue(400);
                }
            }
        });
        busOptions.addNamed("Perform bottom up search", bottomUpSearchEnabled);
        busOptions.addNamed("Perform de novo below m/z", bottomUpSearchOnly);

        parameterPanel.add(new TextHeaderBoxPanel("Bottom Up Search", busOptions));

        card.add(parameterPanel);
        card.add(Box.createRigidArea(new Dimension(0, GuiUtils.SMALL_GAP)));  // gap with element filter

        return card;
    }

    private void setBottomUpSearchMass(JCheckBox bottomUpSearchEnabled) {
        if (bottomUpSearchEnabled.isSelected()) {
            parameterBindings.put("FormulaSearchSettings.enableBottomUpFromMass", () -> "0");
        } else {
            parameterBindings.put("FormulaSearchSettings.enableBottomUpFromMass", () -> String.valueOf(Double.POSITIVE_INFINITY));
        }
    }

//    private JPanel createDeNovoStrategyCard() {
//        enableDeNovo();
//        disableBottomUpSearch();
//
//        return elementFilterPanel;
//    }

    private void disableBottomUpSearch() {
        parameterBindings.put("FormulaSearchSettings.enableBottomUpFromMass", () -> String.valueOf(Double.POSITIVE_INFINITY));
    }

    private void enableDeNovo() {
        parameterBindings.put("FormulaSearchSettings.disableDeNovoAboveMass", () -> String.valueOf(Double.POSITIVE_INFINITY));
    }

    private void disableDeNovo() {
        parameterBindings.put("FormulaSearchSettings.disableDeNovoAboveMass", () -> "0");
    }

    private JPanel createDatabaseStrategyParameters() {
        JPanel card = new JPanel();
        card.setLayout(new BoxLayout(card, BoxLayout.PAGE_AXIS));

        searchDBList = createDatabasePanel();
        searchDBList.setBorder(BorderFactory.createEmptyBorder(0, GuiUtils.LARGE_GAP, 0, 0));

        card.add(searchDBList);
        card.add(Box.createRigidArea(new Dimension(0, GuiUtils.SMALL_GAP)));  // gap with element filter

        disableBottomUpSearch();
        disableDeNovo();

        return card;
    }

    private JCheckboxListPanel<CustomDataSources.Source> createDatabasePanel() {
        if (this.searchDBList != null) return this.searchDBList;
        // configure database to search list
        searchDBList = new JCheckboxListPanel<>(new DBSelectionList(), "Use DB formulas only");
        GuiUtils.assignParameterToolTip(searchDBList.checkBoxList, "FormulaSearchDB");
        parameterBindings.put("FormulaSearchDB", () -> String.join(",", getFormulaSearchDBStrings()));
        PropertyManager.DEFAULTS.createInstanceWithDefaults(FormulaSearchDB.class).searchDBs
                .forEach(s -> searchDBList.checkBoxList.check(CustomDataSources.getSourceFromName(s.name())));
        return searchDBList;
    }

    private JPanel createElementFilterPanel() {
        List<Element> possibleDetectableElements = new ArrayList<>(ApplicationCore.SIRIUS_PROVIDER.sirius().getMs1Preprocessor().getSetOfPredictableElements());
        final FormulaSettings formulaSettings = PropertyManager.DEFAULTS.createInstanceWithDefaults(FormulaSettings.class);

        JPanel panel = applyDefaultLayout(new JPanel());
        JPanel flowPanel = new JPanel();
        final TwoColumnPanel filterFields = new TwoColumnPanel();

        addDefaultStrategyElementFilterSettings(filterFields);
        addDatabaseStrategyElementFilterSettings(filterFields);

        JTextField enforcedTextBox = makeParameterTextField("FormulaSettings.enforced", formulaSettings.getEnforcedAlphabet().toString(), 20);
        final JTextField detectableTextBox;
        filterFields.addNamed("Formula constraints", enforcedTextBox);
        if (isBatchDialog) {
            detectableTextBox = makeParameterTextField("FormulaSettings.detectable", 20);
            filterFields.addNamed("Autodetect", detectableTextBox);
        } else {
            detectableTextBox = null;
        }

        JPanel buttonPanel = new JPanel();

        JButton buttonEdit = new JButton("Customize Filter");
        buttonEdit.setAlignmentX(CENTER_ALIGNMENT);
        buttonPanel.add(buttonEdit);
        filterFields.add((Component) null, buttonPanel);


        buttonEdit.addActionListener(e -> {
            FormulaConstraints currentConstraints = FormulaConstraints.fromString(enforcedTextBox.getText());
            Collection<Element> currentAuto = null;
            if (isBatchDialog) {
                try {
                    currentAuto = ChemicalAlphabet.fromString(detectableTextBox.getText()).toSet();
                } catch (UnknownElementException ex) {
                    currentAuto = possibleDetectableElements;
                }
            }
            ElementSelectionDialog dialog = new ElementSelectionDialog(owner, "Filter Elements", isBatchDialog ? possibleDetectableElements : null, currentAuto, currentConstraints);
            if (dialog.isSuccess()) {
                enforcedTextBox.setText(dialog.getConstraints().toString());
                if (isBatchDialog) {
                    detectableTextBox.setText(join(dialog.getAutoDetect()));
                }
            }
        });

        if (!isBatchDialog) {
            JButton buttonAutodetect = new JButton("Auto detect");
            buttonAutodetect.setAlignmentX(CENTER_ALIGNMENT);
            buttonAutodetect.setToolTipText("Auto detectable element are: " + join(possibleDetectableElements));
            buttonAutodetect.addActionListener(e -> detectElements(enforcedTextBox));
            buttonPanel.add(buttonAutodetect);
        }

        flowPanel.add(filterFields);
        panel.add(new TextHeaderBoxPanel("Element Filter", flowPanel));
        return panel;
    }

    private void addDefaultStrategyElementFilterSettings(TwoColumnPanel filterFields) {
        JComboBox<ElementAlphabetStrategy> elementAlphabetStrategySelector = new JComboBox<>(); //todo NewWorflow: implement this feature in sirius-libs
        List<ElementAlphabetStrategy> settingsElements = java.util.List.copyOf(EnumSet.allOf(ElementAlphabetStrategy.class));
        settingsElements.forEach(elementAlphabetStrategySelector::addItem);
        elementAlphabetStrategySelector.setSelectedItem(ElementAlphabetStrategy.DE_NOVO_ONLY);
        parameterBindings.put("FormulaSearchSettings.applyFormulaContraintsToBottomUp", () -> Boolean.toString(elementAlphabetStrategySelector.getSelectedItem() == ElementAlphabetStrategy.BOTH));

        JLabel label = new JLabel("Apply element filter to");
        filterFields.add(label, elementAlphabetStrategySelector);

        strategyComponents.get(Strategy.DEFAULT).add(label);
        strategyComponents.get(Strategy.DEFAULT).add(elementAlphabetStrategySelector);
    }

    private void addDatabaseStrategyElementFilterSettings(TwoColumnPanel filterFields) {
        JCheckBox useElementFilter = new JCheckBox(); //todo NewWorkflow: implement this feature. This makes the organics filter obsolete. Maybe dont use the checkbox but always select the organics. Make new Element panel popup
        useElementFilter.setSelected(false);
        parameterBindings.put("FormulaSearchSettings.applyFormulaContraintsToCandidateLists", () -> Boolean.toString(useElementFilter.isSelected()));

        JLabel label = new JLabel("Use element filter");
        filterFields.add(label, useElementFilter);

        strategyComponents.get(Strategy.DATABASE).add(label);
        strategyComponents.get(Strategy.DATABASE).add(useElementFilter);

//        elementFilterPanel.setVisible(useElementFilter.isSelected());
//
//        useElementFilter.addActionListener(e -> {
//            elementFilterPanel.setVisible(useElementFilter.isSelected());
//            elementFilterPanel.setEnabled(useElementFilter.isSelected()); //todo ElementFilter: this is not the proper way. Buttons still disabled.
//        });
    }

    private String join(List<?> objects) {
        return objects.stream().map(Object::toString).collect(Collectors.joining(","));
    }

//    protected ElementsPanel createElementPanel(boolean multi) {
//        final FormulaSettings formulaSettings = PropertyManager.DEFAULTS.createInstanceWithDefaults(FormulaSettings.class);
//        List<Element> possDetectableElements = new ArrayList<>(ApplicationCore.SIRIUS_PROVIDER.sirius().getMs1Preprocessor().getSetOfPredictableElements());
//
//        final JButton elementAutoDetect;
//        if (multi) {
//            elementPanel = new ElementsPanel(owner, 4, possDetectableElements, formulaSettings.getAutoDetectionElements(), formulaSettings.getEnforcedAlphabet());
//        } else {
//            /////////////Solo Element//////////////////////
//            elementPanel = new ElementsPanel(owner, 4, formulaSettings.getEnforcedAlphabet());
//            elementAutoDetect = new JButton("Auto detect");
//            elementAutoDetect.setToolTipText("Auto detectable element are: "
//                    + possDetectableElements.stream().map(Element::toString).collect(Collectors.joining(",")));
//            elementAutoDetect.addActionListener(e -> detectElements());
//            elementAutoDetect.setEnabled(true);
//            elementPanel.lowerPanel.add(elementAutoDetect);
//        }
//
//        elementPanel.setBorder(BorderFactory.createEmptyBorder(0, GuiUtils.LARGE_GAP, 0, 0));
//
//        // configure Element panel
//        parameterBindings.put("FormulaSettings.enforced", () -> {
//            return elementPanel.getElementConstraints().toString();
//        });
//        parameterBindings.put("FormulaSettings.detectable", () -> {
//            final List<Element> elementsToAutoDetect = elementPanel.individualAutoDetect ? elementPanel.getElementsToAutoDetect() : Collections.emptyList();
//            return (elementsToAutoDetect.isEmpty() ? "," :
//                    elementsToAutoDetect.stream().map(Element::toString).collect(Collectors.joining(",")));
//        });
//
//        return elementPanel;
//    }

    protected void detectElements(JTextField formulaConstraintsTextBox) {
//        String notWorkingMessage = "Element detection requires MS1 spectrum with isotope pattern.";
//        InstanceBean ec = ecs.get(0);
//        if (!ec.getMs1Spectra().isEmpty() || ec.getMergedMs1Spectrum() != null) {
//            Jobs.runInBackgroundAndLoad(owner, "Detecting Elements...", () -> {
//                final Ms1Preprocessor pp = ApplicationCore.SIRIUS_PROVIDER.sirius().getMs1Preprocessor();
//                ProcessedInput pi = pp.preprocess(new MutableMs2Experiment(ec.getExperiment(), false));
//
//                pi.getAnnotation(FormulaConstraints.class).
//                        ifPresentOrElse(c -> {
//                                    for (Element element : c.getChemicalAlphabet()) {
//                                        if (c.getUpperbound(element) <= 0) {
//                                            c.setLowerbound(element, 0);
//                                            c.setUpperbound(element, 0);
//                                        }
//                                    }
//                                    elementPanel.setSelectedElements(c);
//                                },
//                                () -> new ExceptionDialog(owner, notWorkingMessage)
//                        );
//            }).getResult();
//        } else {
//            new ExceptionDialog(owner, notWorkingMessage);
//        }
    }


    public List<CustomDataSources.Source> getFormulaSearchDBs() {
        return searchDBList.checkBoxList.getCheckedItems();
    }

    public List<String> getFormulaSearchDBStrings() {
        return getFormulaSearchDBs().stream().map(CustomDataSources.Source::id).filter(Objects::nonNull).collect(Collectors.toList());
    }
}
