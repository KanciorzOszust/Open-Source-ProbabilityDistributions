public class MultivariateFishersNoncentralHypergeometric extends MathLibrary{
    int c;
    int[] m;
    int N;
    double n;
    double[] omega;

    MultivariateFishersNoncentralHypergeometric(int[] m, double n, double[] omega) {
        int N = 0;
        for (int i = 0; i < m.length; i++) {
            N += m[i];
        }
        if (n < 0 || n >= N) throw new IllegalArgumentException("0 <= n < sum of m");
        if (omega.length != m.length) throw new IllegalArgumentException("omega.length = m.length");
        for (int i = 0; i < omega.length; i++) {
            if (omega[i] <= 0) throw new IllegalArgumentException("omega[i] > 0");
        }
        this.c = m.length;
        this.m = m;
        this.N = N;
        this.n = n;
        this.omega = omega;
    }

    private double P0support(int y) {
        double value = 1;
        for (int i = 0; i < c; i++) {
            value *= binomial(m[i], y) * power(omega[i], y);
        }
        return value;
    }

    private double P0(int[] S) {
        double value = 0;
        for (int i = 0; i < S.length; i++) {
            value += P0support(S[i]);
        }
        return value;
    }

    double PMF(int[] S) {
        if (S.length != c) throw new IllegalArgumentException("S.length = m.length");
        int testValue = 0;
        for (int i = 0; i < c; i++) {
            testValue += S[i];
        }
        if (testValue != n) throw new IllegalArgumentException("sum of S = n");
        double value = 1;
        for (int i = 0; i < c; i++) {
            value *= binomial(m[i], S[i]) * power(omega[i], S[i]);
        }
        return (1 / P0(S)) * value;
    }
}
