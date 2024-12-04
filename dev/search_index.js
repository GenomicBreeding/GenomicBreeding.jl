var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GenomicBreeding","category":"page"},{"location":"#GenomicBreeding","page":"Home","title":"GenomicBreeding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GenomicBreeding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GenomicBreeding]","category":"page"},{"location":"#GenomicBreeding.Genomes","page":"Home","title":"GenomicBreeding.Genomes","text":"Genomes struct containing unique entries and loci where allele frequencies can have missing values\n\nFields\n\nentries: names of the n entries or samples\nloci: names of the p loci including the chromsome or scaffold name, position, all alleles, and current allele separated by tabs\nallele_frequencies: n x p matrix of allele frequencies between 0 and 1 which can have missing values\nmask: n x p matrix of boolean mask for selective analyses and slicing\n\nExamples\n\njulia> g = Genomes([\"entry_1\", \"entry_2\"], [\"chr1_12345_A|T_A\", \"chr2_678910_C|D_D\"], [0.5 0.25; 0.9 missing], [true true; true false])\nGenomes([\"entry_1\", \"entry_2\"], [\"chr1_12345_A|T_A\", \"chr2_678910_C|D_D\"], Union{Missing, Real}[0.5 0.25; 0.9 missing], Bool[1 1; 1 0])\n\njulia> fieldnames(Genomes)\n(:entries, :loci, :allele_frequencies, :mask)\n\n\n\n\n\n","category":"type"},{"location":"#GenomicBreeding.check-Tuple{Genomes}","page":"Home","title":"GenomicBreeding.check","text":"Check dimension compatibility of the fields of the Genomes struct\n\nExamples\n\njulia> g = Genomes([\"entry_1\", \"entry_2\"], [\"chr1_12345_A|T_A\", \"chr2_678910_C|D_D\"], [0.5 0.25; 0.9 missing], [true true; true false]);\n\njulia> check(g)\ntrue\n\njulia> g.entries = [\"beaking_change\"];\n\njulia> check(g)\nfalse\n\n\n\n\n\n","category":"method"}]
}
