public class Dirichlet extends MathLibrary{
    int K;
    double[] alpha;
    double alpha0;

    Dirichlet(double[] alpha) {
        if (alpha.length < 2) throw new IllegalArgumentException("alpha.length >= 2");
        for (int i = 0; i < alpha.length; i++) {
            if (alpha[i] <= 0) throw new IllegalArgumentException("alpha[i] > 0");
        }
        this.alpha = alpha;
        this.K = alpha.length;
        for (int i = 0; i < alpha.length; i++) {
            this.alpha0 += alpha[i];
        }
    }

    private double B() {
        double value = 1;
        for (int i = 0; i < K; i++) {
            value *= gamma(alpha[i]);
        }
        return value / gamma(alpha0);
    }

    double PDF(double[] x) {
        if (x.length < K) throw new IllegalArgumentException("x.lenth = alpha.length");
        double testValue = 0;
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0 || x[i] > 1) throw new IllegalArgumentException("0 <= x[i] <= 1");
            testValue += x[i];
        }
        if (testValue != 0) throw new IllegalArgumentException("sum of x[i] = 1");

        double value = 1;
        for (int i = 0; i < K; i++) {
            value *= doublePower(x[i], alpha[i] - 1);
        }
        return value / B();
    }

    double Mean(int i) {
        return alpha[i] / alpha0;
    }

    double Mode(int i) {
        if (alpha[i] > 0) {
            return (alpha[i] - 1) / (alpha0 - K);
        } else {
            return Double.NaN;
        }
    }

    private double tildeAlpha(int i) {
        return alpha[i] / alpha0;
    }

    double Variance(int i) {
        return (tildeAlpha(i) * (1 - tildeAlpha(i))) / (alpha0 + 1);
    }

    double Entropy() {
        double value = 0;
        for (int i = 0; i < K; i++) {
            value += (alpha[i] - 1) * digamma(alpha[i]);
        }
        return ln(B()) + ((alpha0 - K) * digamma(alpha0)) - value;
    }

    double MethodOfMoments(int i, int j) {
        double part = (Mean(j) * (1 - Mean(j))) / Variance(j) - 1;
        return Mean(i) * part;
    }
}
