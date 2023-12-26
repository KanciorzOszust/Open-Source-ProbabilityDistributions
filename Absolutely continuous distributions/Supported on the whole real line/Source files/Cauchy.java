public class Cauchy extends MathLibrary{
    double x0;
    double gamma;

    Cauchy(double x0, double gamma) {
        if (gamma <= 0) throw new IllegalArgumentException("gamma > 0");
        this.x0 = x0;
        this.gamma = gamma;
    }

    double PDF(double x) {
        return 1 / (PI * gamma * (1 + power((x - x0) / gamma, 2)));
    }

    double CDF(double x) {
        return (1 / PI) * arcTangent((x - x0) / gamma) + 0.5;
    }

    double Quantile(double p) {
        return x0 + (gamma * tangent(PI * (p - 0.5)));
    }

    double Median() {
        return x0;
    }

    double Mode() {
        return x0;
    }

    double MAD() {
        return gamma;
    }

    double Entropy() {
        return log(4 * PI * gamma, 10);
    }

    double FisherInformation() {
        return 1 / (2 * power(gamma, 2));
    }
}
