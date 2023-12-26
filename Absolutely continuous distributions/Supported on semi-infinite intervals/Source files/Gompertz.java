public class Gompertz extends MathLibrary{
    double eta;
    double b;

    Gompertz(double eta, double b) {
        if (eta <= 0 || b <= 0) throw new IllegalArgumentException("eta, b > 0");
        this.eta = eta;
        this.b = b;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return b * eta * euler(eta + (b * x) - (eta * euler(b * x)));
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - euler(-1 * eta * (euler(b * x) - 1));
    }

    double Mean() {
        return (1 / b) * euler(eta) * exponentialIntergral(-1 * eta);
    }

    double Median() {
        return (1 / b) * ln((1 / eta) * ln(0.5) + 1);
    }

    double Mode() {
        if (eta >= 1) {
            return 0;
        } else {
            return (1 / b) * ln(1 / eta);
        }
    }
}
