package picard.vcf;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * A class to store the various classifications for:
 * 1. a truth genotype versus a reference
 * 2. a call genotype versus a truth, relative to a reference
 *
 * An example use of this class is to have one instance per following use case:
 * - SNPs
 * - indels
 * - filtered variant (truth or call)
 * - filtered genotype (truth or call)
 * - low GQ (call)
 * - low DP (call)
 * - No call (truth or call)
 * - No variant (truth or call) *
 *
 * @author nhomer
 */
public class GenotypeConcordanceStates {
    /**
     * These states represent the relationship between a truth genotype and the reference sequence.
     */
    public enum TruthState {
        HOM_REF, // ref/ref
        HET_REF_VAR1, // ref/var1 (var1!=ref)
        HET_VAR1_VAR2, // var1/var2 (var1!=var2, var1!=ref, var2!=ref)
        HOM_VAR1, // var1/var1 (var1!=ref)
        NO_CALL,
        LOW_GQ,
        LOW_DP,
        FILTERED,
        IS_MIXED
    }

    /**
     * These states represent the relationship between the call genotype and the truth genotype relative to
     * a reference sequence.
     */
    enum CallState {
        HOM_REF, // ref/ref, valid for all TruthStates
        HET_REF_VAR1, // ref/var1, valid for all TruthStates
        HET_REF_VAR2, // ref/var2, valid only for TruthStates: HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1
        HET_REF_VAR3, // ref/var3, valid only for TruthStates: HET_VAR1_VAR2
        HET_VAR1_VAR2, // var1/var2, valid for all TruthStates
        HET_VAR1_VAR3, // var1/var3, valid only for TruthStates: HET_VAR1_VAR2
        HET_VAR3_VAR4, // var3/var4, valid only for TruthStates: HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1
        HOM_VAR1, // var1/var1, valid for all TruthStates
        HOM_VAR2, // var2/var2, valid only for TruthStates: HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1
        HOM_VAR3, // var3/var3, valid only for TruthStates: HET_VAR1_VAR2
        NO_CALL,
        LOW_GQ,
        LOW_DP,
        FILTERED,
        IS_MIXED
    }

    static final Set<TruthState> allHetTruthStates = new HashSet<TruthState>(Arrays.asList(TruthState.HET_REF_VAR1, TruthState.HET_VAR1_VAR2));
    static final Set<TruthState> allHomVarTruthStates = new HashSet<TruthState>(Arrays.asList(TruthState.HOM_VAR1));
    static final Set<TruthState> allVarTruthStates = new HashSet<TruthState>();

    static final Set<CallState> allHetCallStates = new HashSet<CallState>(Arrays.asList(CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3,
                                                                                        CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4));
    static final Set<CallState> allHomVarCallStates = new HashSet<CallState>(Arrays.asList(CallState.HOM_VAR1, CallState.HOM_VAR2, CallState.HOM_VAR3));
    static final Set<CallState> allVarCallStates = new HashSet<CallState>();

    static {
        allVarTruthStates.addAll(allHetTruthStates);
        allVarTruthStates.addAll(allHomVarTruthStates);

        allVarCallStates.addAll(allHetCallStates);
        allVarCallStates.addAll(allHomVarCallStates);
    }






    /**
     * A specific state for a 2x2 contingency table.
     */
    enum ContingencyState {
        TP,
        FP,
        TN,
        FN,
        NA
    }

    /**
     * A minute class to store the truth and call state respectively.
     */
    static class TruthAndCallStates implements Comparable<TruthAndCallStates>{
        public final TruthState truthState;
        public final CallState callState;

        public TruthAndCallStates(final TruthState truthState, final CallState callState) {
            this.truthState = truthState;
            this.callState = callState;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            return compareTo((TruthAndCallStates) o) == 0;
        }

        @Override
        public int hashCode() {
            int result = truthState.hashCode();
            result = 31 * result + callState.hashCode();
            return result;
        }

        @Override
        public int compareTo(final TruthAndCallStates that) {
            int result = this.truthState.compareTo(that.truthState);
            if (result == 0) result = this.callState.compareTo(that.callState);
            return result;
        }
    }
}
