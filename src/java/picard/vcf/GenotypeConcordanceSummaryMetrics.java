package picard.vcf;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Class that holds summary metrics about Genotype Concordance
 *
 * @author George Grant
 */
public class GenotypeConcordanceSummaryMetrics extends MetricBase {
    GenotypeConcordanceSummaryMetrics(final String eventType, final GenotypeConcordanceCounts concordanceCounts,
                                      final String truthSampleName, final String callSampleName) {
        this.EVENT_TYPE = eventType;
        this.TRUTH_SAMPLE_NAME = truthSampleName;
        this.CALL_SAMPLE_NAME = callSampleName;

        final GenotypeConcordanceScheme scheme = new GenotypeConcordanceScheme();
        concordanceCounts.validateCountsAgainstScheme(scheme);

        this.HET_SENSITIVITY = concordanceCounts.getSensitivity(scheme, GenotypeConcordanceStates.allHetTruthStates);
        this.HET_PPV = concordanceCounts.getPPV(scheme, GenotypeConcordanceStates.allHetCallStates);

        this.HOMVAR_SENSITIVITY = concordanceCounts.getSensitivity(scheme, GenotypeConcordanceStates.allHomVarTruthStates);
        this.HOMVAR_PPV = concordanceCounts.getPPV(scheme, GenotypeConcordanceStates.allHomVarCallStates);

        this.VAR_SENSITIVITY = concordanceCounts.getSensitivity(scheme, GenotypeConcordanceStates.allVarTruthStates);
        this.VAR_PPV = concordanceCounts.getPPV(scheme, GenotypeConcordanceStates.allVarCallStates);
    }

    /** The type of event SNP/Indel */
    public String EVENT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE_NAME;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE_NAME;

    /** The het sensitivity */
    public double HET_SENSITIVITY;

    /** The het ppv (positive predictive value) */
    public double HET_PPV;

    /** The homVar sensitivity */
    public double HOMVAR_SENSITIVITY;

    /** The homVar ppv (positive predictive value) */
    public double HOMVAR_PPV;

    /** The var sensitivity */
    public double VAR_SENSITIVITY;

    /** The var ppv (positive predictive value) */
    public double VAR_PPV;
}
