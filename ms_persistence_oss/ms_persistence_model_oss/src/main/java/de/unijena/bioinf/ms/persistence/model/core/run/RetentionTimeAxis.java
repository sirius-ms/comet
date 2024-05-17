package de.unijena.bioinf.ms.persistence.model.core.run;

import jakarta.persistence.Id;
import lombok.*;
import lombok.experimental.SuperBuilder;

/**
 * Each LC/MS run maps an index to a scan number and a retention time
 */
@Getter
@Setter
@SuperBuilder
@AllArgsConstructor
@NoArgsConstructor
public class RetentionTimeAxis {

    @Id
    private long runId;
    private int[] scanIndizes;
    private int[] scanNumbers;
    private String[] scanIdentifiers;
    private double[] retentionTimes;
    /**
     * Estimated noise level per scan to compute a local SNR value.
     */
    private float[] noiseLevelPerScan;

}
