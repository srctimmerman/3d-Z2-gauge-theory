using Plots
using DelimitedFiles
using LaTeXStrings

#loads data for various lattice sizes
four = readdlm("4.txt")
six = readdlm("6.txt")
eight = readdlm("8.txt")
ten = readdlm("10.txt")

#plots the specific heat as a function of temperature for each of the lattice
# sizes
plot([four[:,1], six[:,1], eight[:,1], ten[:,1]], [four[:,3], six[:,3], eight[:,3], ten[:,3]],
    markershape= :x, markersize = 2, label = [L"L=4" L"L=6" L"L=8" L"L=10" L"L=12"],
    xticks = ([1.00, 1.3133, 1.50, 2.00], [L"1.00", L"1.3133", L"1.50", L"2.00"]),
    yticks = ([0.0, 0.5, 1.0, 1.5, 2.0], [L"0.0", L"0.5", L"1.0", L"1.5", L"2.0"]);
    palette = palette([:blue, :orange], 4))

plot!([1/0.7614125], seriestype="vline", seriescolor = "gray", linewidth = 2, linestyle = :dot, label = false)
