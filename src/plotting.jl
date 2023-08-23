function plot_querycover(
    query_cover::Vector{Float64},
    fileout::AbstractString
)
    categories = [(0.0, 20.99), (21.0, 40.99), (41.0, 60.99), (61.0, 80.99), (81.0, 100.0)]
    labels = ["0.01-20%", "21-40%", "41-60%", "61-80%", "81-100%"]
    bincount = map(category ->
        count(x -> first(category) <= x <= last(category), query_cover),
        categories)
    toptitle = "Query coverage in markers"
    println(toptitle, ":\n", UnicodePlots.barplot(labels, bincount))
    Plots.bar(labels, bincount, color = "green", title = toptitle, bar_width = 0.3,
        legend = false)
    Plots.savefig(fileout)
    return fileout
end

function plot_percentident(
    identities::Vector{Float64},
    subjectcount::Int,
    fileout::AbstractString
)
    querycount = length(identities)
    @assert subjectcount > 0
    @assert subjectcount >= querycount
    categories = [(0.0, 30.99), (31.0, 50.99), (51.0, 70.99), (71.0, 90.99), (91.0, 100.0)]
    values = map(category ->
        count(x -> first(category) <= x <= last(category), identities),
        categories)
    labels = ["missing", "≤ 30%", "> 31-50%", "> 51-70%", "> 71-90%", "> 91-100%"]
    missingmarkers = subjectcount - querycount
    pushfirst!(values, missingmarkers)
    toptitle = "Percent identity in markers"
    println(toptitle, ":\n", UnicodePlots.barplot(labels, values))
    Plots.pie(labels, values, palette = :Greens_8, title = toptitle, linecolor = "black")
    Plots.savefig(fileout)
    return fileout
end

