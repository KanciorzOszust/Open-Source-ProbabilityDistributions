public class GammaGompertz extends MathLibrary{
    double b;
    double s;
    double beta;

    GammaGompertz(double b, double s, double beta) {
        if (b <= 0 || s<= 0 || beta <= 0) throw new IllegalArgumentException("b, s, beta > 0");
        this.b = b;
        this.s = s;
        this.beta = beta;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = b * s * euler(b * x) * doublePower(beta, s);
        double part2 = doublePower(beta - 1 + euler(b * x), s + 1);
        return part1 / part2;
    }

    double CDF(double x) {
        if (beta == 1) {
            return 1 - euler(-1 * b * s * x);
        } else {
            double part = doublePower(beta - 1 + euler(b * s), s);
            return (1 - doublePower(beta, s)) / part;
        }
    }

    double Mean() {
        if (beta == 1) {
            return 1 / (b * s);
        } else if (s == 1) {
            return (1 / b) * (beta / (beta - 1)) * ln(beta);
        } else {
            return (1 / b) * (1 / s) * hypergeometric(s, 1, s + 1, (beta - 1) / beta);
        }
    }

    double Median() {
        return (1 / b) * ln(beta * (doublePower(0.5, -1 / s) - 1) + 1);
    }

    double Mode() {
        if (beta <= s + 1) {
            return 0;
        } else {
            return (1 / b) * ln((1 / s) * (beta - 1));
        }
    }

    double MGF(double t) {
        if (beta == 1) {
            return (s * b) / (t + (s * b));
        } else {
            double part = doublePower(beta, s) * ((s * b) / (t + (s * b)));
            return part * hypergeometric(s + 1, (t / b) + s, (t / b) + s + 1, 1 - beta);
        }
    }
}
