public class FishersZ extends MathLibrary{
    double d1;
    double d2;

    FishersZ(double d1, double d2) {
        if (d1 <= 0 || d2 <= 0) throw new IllegalArgumentException("d1, d2 > 0");
        this.d1 = d1;
        this.d2 = d2;
    }

    double PDF(double x) {
        double part1 = (2 * doublePower(d1, d1 / 2) * doublePower(d2, d2 / 2)) / beta(d1 / 2, d2 / 2);
        double part2 = euler(d1 * x) / doublePower(d1 * euler(2 * x) + d2, (d1 + d2) / 2);
        return part1 * part2;
    }

    int Mode() {
        return 0;
    }
}
