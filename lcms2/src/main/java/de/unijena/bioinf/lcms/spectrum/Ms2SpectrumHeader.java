package de.unijena.bioinf.lcms.spectrum;

import de.unijena.bioinf.ChemistryBase.ms.CollisionEnergy;
import de.unijena.bioinf.ChemistryBase.ms.IsolationWindow;
import lombok.Getter;

import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.Optional;

public class Ms2SpectrumHeader extends Ms1SpectrumHeader implements Serializable {
    @Nullable protected final CollisionEnergy energy;

    @Nullable protected final IsolationWindow isolationWindow;
    @Getter protected final int parentId, parentIndex;
    @Getter  protected final double retentionTime;
    @Getter protected final double precursorMz;
    @Getter protected final double targetedMz;
    @Getter protected final int msLevel;

    public Ms2SpectrumHeader(String sourceId, int index, int polarity, int mslevel, boolean centroided, CollisionEnergy energy, IsolationWindow window, int parentId, int parentIndex, double precursorMz, double targetedMz, double retentionTime) {
        this(-1, index, sourceId, polarity,mslevel,centroided,energy,window,parentId,parentIndex,precursorMz,targetedMz,retentionTime);
    }

    public Optional<CollisionEnergy> getEnergy() {
        return Optional.ofNullable(energy);
    }

    public Optional<IsolationWindow> getIsolationWindow() {
        return Optional.ofNullable(isolationWindow);
    }

    public Ms2SpectrumHeader(int uid, int index, String sourceId, int polarity, int msLevel, boolean centroided,  CollisionEnergy energy, IsolationWindow window, int parentId, int parentIndex, double precursorMz, double targetedMz, double retentionTime) {
        super(uid, index, sourceId, polarity, centroided);
        this.energy = energy;
        this.parentId = parentId;
        this.isolationWindow = window==null || window.isUndefined() ? null : window;
        this.retentionTime = retentionTime;
        this.parentIndex = parentIndex;
        this.msLevel = msLevel;
        if (precursorMz<0 && targetedMz<0) {
            throw new IllegalArgumentException("Neither precursor nor targeted mz is known");
        }
        if (window!=null) {
            this.precursorMz = precursorMz>=0 ? precursorMz : targetedMz + window.getWindowOffset();
            this.targetedMz = targetedMz>=0 ? targetedMz : precursorMz - window.getWindowOffset();
        } else {
            this.precursorMz = precursorMz;
            this.targetedMz = targetedMz;
        }
    }

    public Ms2SpectrumHeader withUid(int uid) {
        return new Ms2SpectrumHeader(uid, scanIndex, sourceId, polarity, msLevel, centroided, energy, isolationWindow, parentId, parentIndex, precursorMz, targetedMz, retentionTime);
    }

    @Override
    public String toString() {
        return "Ms2SpectrumHeader{" +
                "parentId=" + parentId +
                ", retentionTime=" + retentionTime +
                ", precursorMz=" + precursorMz +
                '}';
    }
}
