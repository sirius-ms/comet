package de.unijena.bioinf.ms.gui.fingerid;

import de.unijena.bioinf.chemdb.DataSource;
import de.unijena.bioinf.chemdb.custom.CustomDataSources;
import de.unijena.bioinf.ms.gui.table.ActiveElementChangedListener;
import de.unijena.bioinf.ms.gui.utils.WrapLayout;
import de.unijena.bioinf.projectspace.FormulaResultBean;

import javax.swing.*;
import java.awt.*;
import java.util.List;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;

public class DBFilterPanel extends JPanel implements ActiveElementChangedListener<FingerprintCandidateBean, Set<FormulaResultBean>>, CustomDataSources.DataSourceChangeListener {
    public final static Set<String> BLACK_LIST = Set.of(DataSource.ADDITIONAL.realName, DataSource.ALL.realName, DataSource.ALL_BUT_INSILICO.realName,
            DataSource.PUBCHEMANNOTATIONBIO.realName, DataSource.PUBCHEMANNOTATIONDRUG.realName, DataSource.PUBCHEMANNOTATIONFOOD.realName, DataSource.PUBCHEMANNOTATIONSAFETYANDTOXIC.realName,
            DataSource.SUPERNATURAL.realName
    );

    private final List<FilterChangeListener> listeners = new LinkedList<>();

    protected long bitSet;
    protected List<JCheckBox> checkboxes;
    private final AtomicBoolean isRefreshing = new AtomicBoolean(false);


    public DBFilterPanel(StructureList sourceList) {
        setLayout(new WrapLayout(FlowLayout.LEFT, 5, 1));
        this.checkboxes = new ArrayList<>(CustomDataSources.size());
        for (CustomDataSources.Source source : CustomDataSources.sources()) {
            if (!BLACK_LIST.contains(source.name()))
                checkboxes.add(new JCheckBox(source.name()));
        }
        addBoxes();
        CustomDataSources.addListener(this);
        sourceList.addActiveResultChangedListener(this);
    }

    public void addFilterChangeListener(FilterChangeListener listener) {
        listeners.add(listener);
    }

    public void fireFilterChangeEvent() {
        listeners.forEach(l -> l.fireFilterChanged(bitSet));
    }

    protected void addBoxes() {
        checkboxes.sort(Comparator.comparing(o -> o.getText().toUpperCase()));

        this.bitSet = 0L;
        for (final JCheckBox box : checkboxes) {
            if (box.isSelected())
                this.bitSet |= CustomDataSources.getSourceFromName(box.getText()).flag();
            add(box);
            box.addChangeListener(e -> {
                if (!isRefreshing.get()) {
                    if (box.isSelected())
                        bitSet |= CustomDataSources.getSourceFromName(box.getText()).flag();
                    else
                        bitSet &= ~CustomDataSources.getSourceFromName(box.getText()).flag();
                    fireFilterChangeEvent();
                }
            });
        }
    }

    protected void reset() {
        isRefreshing.set(true);
        bitSet = 0;
        try {
            for (JCheckBox checkbox : checkboxes) {
                checkbox.setSelected(false);
            }
        } finally {
            fireFilterChangeEvent();
            isRefreshing.set(false);
        }
    }

    public boolean toggle() {
        setVisible(!isVisible());
        return isVisible();
    }

    @Override
    public void resultsChanged(Set<FormulaResultBean> datas, FingerprintCandidateBean sre, List<FingerprintCandidateBean> resultElements, ListSelectionModel selections) {
        reset();
    }

    @Override
    public void fireDataSourceChanged(Collection<String> changes) {
        HashSet<String> changed = new HashSet<>(changes);
        isRefreshing.set(true);
        boolean c = false;
        Iterator<JCheckBox> it = checkboxes.iterator();

        while (it.hasNext()) {
            JCheckBox checkbox = it.next();
            if (changed.remove(checkbox.getText())) {
                it.remove();
                c = true;
            }
        }

        for (String name : changed) {
            checkboxes.add(new JCheckBox(name));
            c = true;
        }

        if (c) {
            removeAll();
            addBoxes();
            revalidate();
            repaint();
            fireFilterChangeEvent();
        }


        isRefreshing.set(false);
    }

    public interface FilterChangeListener extends EventListener {
        void fireFilterChanged(long filterSet);
    }
}
