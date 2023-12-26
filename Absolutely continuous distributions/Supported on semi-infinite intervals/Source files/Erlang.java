public class Erlang extends MathLibrary{
    int k;
    double gamma;

    Erlang(int k, double gamma) {
        if (k < 1) throw new IllegalArgumentException("k >= 1");
        if (gamma <= 0) throw new IllegalArgumentException("gamma > 0");
        this.k = k;
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part = power(gamma, k) * power(x, k - 1) * euler(-1 * gamma * x);
        return part / factorial(k - 1);
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return incompleteLowerGamma(k, gamma * x) / factorial(k - 1);
    }

    double Mean() {
        return k / gamma;
    }

    double Mode() {
        return (1 / gamma) * (k - 1);
    }

    double Variance() {
        return k / power(gamma, 2);
    }

    double Skewness() {
        return 2 / root(k, 2);
    }

    double EXkurtosis() {
        return 6.0 / k;
    }

    double Entropy() {
        return (1 - k) * digamma(k) + ln(gamma(k) / gamma) + k;
    }

    double MGF(double t) {
        if (t < gamma) {
            return doublePower(1 - (t / gamma), -1 * k);
        } else {
            return Double.NaN;
        }
    }
}
