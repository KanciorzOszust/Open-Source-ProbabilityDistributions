public class RaisedCosine extends MathLibrary{
    double mu;
    double s;

    RaisedCosine(double mu, double s) {
        if (s <= 0) throw new IllegalArgumentException("s > 0");
        this.mu = mu;
        this.s = s;
    }

    double PDF(double x) {
        if (x < (mu - s) || x > (mu + s)) throw new IllegalArgumentException("mu - s <= x <= mu + a");
        return (1 / (2 * s)) * (1 + cosine((x - mu) / s * PI));
    }

    double CDF(double x) {
        if (x < (mu - s) || x > (mu + s)) throw new IllegalArgumentException("mu - s <= x <= mu + a");
        return 0.5 * (1 + ((x - mu) / s) + ((1 / PI) * sine((x - mu) / s * PI)));
    }

    double Mean() {
        return mu;
    }

    double Median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return power(s, 2) * ((1 / 3.0) - (2 / power(PI, 2)));
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return (6 * (90 - power(PI, 4))) / (5 * power(power(PI, 2) - 6, 2));
    }

    double MGF(double t) {
        double part1 = power(PI, 2) * hyperbolicSine(s * t);
        double part2 = s * t * (power(PI, 2) - (power(s, 2) * power(t, 2)));
        return part1 / part2 * euler(mu * t);
    }
}
