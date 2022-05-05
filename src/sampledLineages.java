import java.util.List;

public class sampledLineages {

    List<Integer> lineages;
    List<Double> times;

    public sampledLineages(List<Integer> lineages, List<Double> times) {

        this.lineages = lineages;
        this.times = times;

    }

    public List<Integer> getSampledLineages() {

        return this.lineages;
    }

    public List<Double> getSampledTimes() {

        return this.times;
    }

}

