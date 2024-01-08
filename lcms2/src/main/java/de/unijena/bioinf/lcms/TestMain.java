package de.unijena.bioinf.lcms;

import de.unijena.bioinf.ChemistryBase.jobs.SiriusJobs;
import de.unijena.bioinf.ChemistryBase.math.Statistics;
import de.unijena.bioinf.jjobs.BasicJJob;
import de.unijena.bioinf.jjobs.JJob;
import de.unijena.bioinf.jjobs.JobManager;
import de.unijena.bioinf.lcms.align2.AlignmentBackbone;
import de.unijena.bioinf.lcms.merge2.MergedTrace;
import de.unijena.bioinf.lcms.trace.ProcessedSample;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.*;
import java.util.*;

/**
 * Aktuelle Vorgehensweise:
 * - erstmal alle Samples und deren Traces extrahieren und in die DB speichern
 * - Traces werden dann nochmal in "Trace Chains" organisiert, wobei eine Chain alle Traces mit ähnlicher Masse und
 *   unterschiedlicher (nicht überlappender!) Retentionszeit enthält. Dieser Schritt ist nicht notwendig und kann
 *   wieder rausgenommen werden. Der einzige Vorteil von diesem Schritt ist eigentlich, dass man die minimale
 *   Massenabweichung bekommt, die zwei klar unterschiedliche Traces haben dürfen.
 *
 * - danach werden die Apexe von jedem Trace gesammelt und aligniert
 *      - Alignment findet stufenweise statt
 *      - zuerst aligniert man Apexe die sehr gut aussehen (=Isotopenpeaks haben und/oder hohe Intensität)
 *      - danach wird rekalibriert und nochmal neu aligniert, diesmal alle Apexe
 *      - für jeden Apex speichern wir das "Rechteck" ab, in dem sein Trace sich befindet. D.h. wir wissen
 *        die m/z und rt range über die der Trace verläuft
 * - die Rekalibrierung dient erstmal nur dem bestimmen der Rekalibrierungsfunktionen für m/z und rt. m/z
 *   Rekalibrierung scheint auf den Testdaten nichts zu bringen, aber wer weiß
 *
 * - im nächsten Schritt gehen wir über alle gemergten Apexe und bestimmen die Vereinigung aller zugehörigen Rechtecke
 * - liegen zwei Rechtecke ineinander or haben geringe Massenabweichung werden sie gemerged
 * - ansonsten werden sie separiert, einmal in m/z Richtung und einmal in rt Richtung. So bekommt man zwei Rechtecke,
 *   eines ist breiter in der Retentionszeit, hat aber geringere Massenabweichung, eins ist breiter in der Massenabweichung,
 *   hat aber geringe RT Zeit
 * - alle Rechtecke sind jetzt disjunkt, wir können also nochmal über alle Samples durchgehen und jedes Rechteck nehmen
 *   und alle Intensitäten darin aufsummieren. Für die "Doppel-Rechtecke" gehen wir über beides drüber (sammeln also Peaks
 *   im engen Retentionszeitfenster mit höherer Massenabweichung ein und dann nochmal die äußeren Peaks mit geringerer
 *   Massenabweichung).
 *
 * - ob die doppelten Rechtecke sinnvol sind? Keine Ahnung, sie erlauben aber jedenfalls dass wir am Ende klar definierte
 *   Regionen samplen können, was wiederum den Vorteil hat, dass wir nie versehentlich zwei Peaks doppelt samplen.
 *
 *
 *
 */
public class TestMain {

    public static void main(String[] args) {
        final de.unijena.bioinf.lcms.trace.ProcessedSample[] samples;
        LCMSProcessing processing = new LCMSProcessing();
        {
            JobManager globalJobManager = SiriusJobs.getGlobalJobManager();
            //File[] files = new File("/home/kaidu/analysis/lcms/diverse_collection/small").listFiles();
            File[] files = new File("/home/kaidu/analysis/lcms/diverse_collection/MSV000080627/").listFiles();
            //File[] files = new File("/home/kaidu/data/raw/polluted_citrus/").listFiles();
            List<BasicJJob<de.unijena.bioinf.lcms.trace.ProcessedSample>> jobs = new ArrayList<>();
            int atmost = Integer.MAX_VALUE;
            for (File f : files) {
                if (--atmost < 0) break;
                if (f.getName().endsWith(".mzML")) {
                    jobs.add(SiriusJobs.getGlobalJobManager().submitJob(new BasicJJob<de.unijena.bioinf.lcms.trace.ProcessedSample>() {
                        @Override
                        protected de.unijena.bioinf.lcms.trace.ProcessedSample compute() throws Exception {
                            ProcessedSample sample = processing.processSample(f);
                            sample.inactive();
                            return sample;
                        }
                    }));
                }
            }
            //samples = jobs.stream().map(JJob::takeResult).toArray(de.unijena.bioinf.lcms.trace.ProcessedSample[]::new);
            for (BasicJJob<ProcessedSample> job : jobs) {
                System.out.println(job.takeResult().getUid());
            }
        }
        try {
            AlignmentBackbone bac = processing.align();
            ProcessedSample merged = processing.merge(bac);
            DoubleArrayList avgAl = new DoubleArrayList();
            for (MergedTrace r : merged.getTraceStorage().getMergeStorage()) {
                System.out.println(r);
                avgAl.add(r.getSampleIds().size());
            }
            System.out.println("AVERAGE = " + avgAl.doubleStream().sum()/avgAl.size());
            System.out.println("Good Traces = " + avgAl.doubleStream().filter(x->x>=5).sum());
            System.out.println("-----------------------------");
            processing.exportFeaturesToFiles(merged, bac);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }
}
