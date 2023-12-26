public class ShiftedGompertz extends MathLibrary{
    double b;
    double eta;

    ShiftedGompertz(double b, double eta) {
        if (b < 0 || eta < 0) throw new IllegalArgumentException("b, eta >= 0");
        this.b = b;
        this.eta = eta;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = b * euler(-1 * b * x) * euler(-1 * eta * euler(-1 * b * x));
        double part2 = 1 + (eta * (1 - euler(-1 * b * x)));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (1 - euler(-1 * b * x)) * euler(-1 * eta * euler(-1 * b * x));
    }

    double Mode() {
        if (eta > 0.5) {
            double part = (3 + eta - root(power(eta, 2) + (2 * eta) + 5, 0)) / (2 * eta);
            return (-11 / power(b, 2)) * ln(part);
        } else {
            return 0;
        }
    }
}
