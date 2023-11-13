package de.unijena.bioinf.lcms.spectrum;

import java.io.Serializable;

public class Ms1SpectrumHeader implements Serializable {

    protected final int uid;
    protected final byte polarity;

    protected final boolean centroided;

    public Ms1SpectrumHeader(int polarity, boolean centroided) {
        this(-1, polarity, centroided);
    }

    public Ms1SpectrumHeader(int uid, int polarity, boolean centroided) {
        this.uid = uid;
        this.polarity = (byte)polarity;
        this.centroided = centroided;
    }

    public int getUid() {
        return uid;
    }

    public Ms1SpectrumHeader withUid(int uid) {
        return new Ms1SpectrumHeader(uid, polarity, centroided);
    }
}
