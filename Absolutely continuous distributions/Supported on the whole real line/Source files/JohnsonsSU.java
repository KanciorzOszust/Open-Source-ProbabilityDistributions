public class JohnsonsSU extends MathLibrary{
    double gamma;
    double xi;
    double delta;
    double lambda;

    JohnsonsSU(double gamma, double xi, double delta, double lambda) {
        if (gamma <= 0 || xi <= 0 || delta <= 0 || lambda <= 0) throw new IllegalArgumentException("arguments > 0");
        this.gamma = gamma;
        this.xi = xi;
        this.delta = delta;
        this.lambda = lambda;
    }

    double PDF(double x) {
        double part1 = delta / (lambda * root(2 * PI, 2));
        double part2 = 1 / root(1 + power((x - xi) / lambda, 2), 2);
        double part3 = gamma + (delta * hyperbolicArcSine((x - xi) / lambda));
        return part1 * part2 * euler(-0.5 * power(part3, 2));
    }

    private double PHI(double x) {
        return 0.5 * (1 + errorFunction((x - xi) / lambda));
    }

    double CDF(double x) {
        return PHI(gamma + (delta * hyperbolicArcSine((x - xi) / lambda)));
    }

    double Mean() {
        return xi - (lambda * euler(power(delta, -2) / 2) * hyperbolicSine(gamma / delta));
    }

    double Median() {
        return xi + (lambda * hyperbolicSine(-1 * gamma / delta));
    }

    double Variance() {
        double part1 = (power(lambda, 2) / 2) * (euler(power(delta, -2)) - 1);
        double part2 = euler(power(delta, -2)) * hyperbolicCosine(2 * gamma / delta) + 1;
        return part1 * part2;
    }

    double Skewness() {
        double part1 = power(lambda, 3) * euler(1 / delta) * power(euler(power(delta, -1)) - 1, 2);
        double part2 = euler(power(delta, -2)) * (euler(power(delta, -2)) + 2) * hyperbolicSine(3 * gamma / delta);
        return (part1 * part2 + (3 * hyperbolicSine(2 * gamma / delta))) / (4 * doublePower(Variance(), 1.5));
    }

    double EXkurtosis() {
        double K1part = euler(power(delta, 2)) + (2 * euler(delta)) + (3 * euler(1)) - 3;
        double K1 = euler(1) * K1part * hyperbolicCosine(4 * gamma / delta);
        double K2part = 4 * euler(1) * (euler(power(delta, -2)) + 2);
        double K2 = K2part * hyperbolicCosine(3 * gamma / delta);
        double K3 = 3 * (2 * euler(power(delta, -1)) + 1);
        double part = power(lambda, 4) * power((euler(power(delta, -1)) - 1), 2);
        return (part * (K1 + K2 + K3)) / (8 * power(Variance(), 2));
    }
}
