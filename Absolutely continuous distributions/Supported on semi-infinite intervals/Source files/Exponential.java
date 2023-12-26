public class Exponential extends MathLibrary{
    double gamma;

    Exponential(double gamma) {
        if (gamma <= 0) throw new IllegalArgumentException("gamma > 0");
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return gamma * euler(-1 * gamma * x);
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - euler(-1 * gamma * x);
    }

    double Quantile(double p) {
        return -1 * ln(1 - p) / gamma;
    }

    double Mean() {
        return 1 / gamma;
    }

    double Median() {
        return ln(2) / gamma;
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        return 1 / power(gamma, 2);
    }

    int Skewness() {
        return 2;
    }

    int EXkurtosis() {
        return 6;
    }

    double Entropy() {
        return 1 - ln(gamma);
    }

    double MGF(double t) {
        if (t < gamma) return gamma / (gamma - t);
        else return Double.NaN;
    }

    double FisherInformation() {
        return Variance();
    }

    double KullbackLeiblerDivergence(double gamma0) {
        return ln(gamma0 / gamma) + (gamma / gamma0) - 1;
    }

    double ExpectedShortfall(double p) {
        return -1 * (ln(1 - p) + 1) / gamma;
    }
}
