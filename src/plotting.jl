import Plots
import UnicodePlots

function plot_querycover(
    query_cover::Vector{Float64},
    fileout::AbstractString
)
    categories = [(0.0, 20.0), (21.0, 40.0), (41.0, 60.0), (61.0, 80.0), (81.0, 100.0)]
    labels = ["0.01-20%", "21-40%", "41-60%", "61-80%", "81-100%"]
    bincount = map(category ->
            count(x -> first(category) <= x <= last(category),
                query_cover),
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
    aboveninety = length(eachrow(filter(identity -> identity > 90.0, identities)))
    aboveseventy = length(eachrow(filter(identity -> identity > 70.0, identities)))
    abovefifty = length(eachrow(filter(identity -> identity > 50.0, identities)))
    abovethirty = length(eachrow(filter(identity -> identity > 30.0, identities)))
    belowthirty = querycount - abovefifty
    labels = ["missing", "≤ 30%", "> 30%", "> 50%", "> 70%", "> 90%"]

    missing = round(((subjectcount - querycount) / subjectcount) * 100.0, digits = 2)
    percent_aboveninety = round((aboveninety / subjectcount) * 100.0, digits = 2)
    percent_aboveseventy = round((aboveseventy / subjectcount) * 100.0, digits = 2)
    percent_abovefifty = round((abovefifty / subjectcount) * 100.0, digits = 2)
    percent_abovethirty = round((abovethirty / subjectcount) * 100.0, digits = 2)
    percent_belowthirty = round((belowthirty / subjectcount) * 100.0, digits = 2)
    values = [missing, percent_belowthirty, percent_abovethirty, percent_abovefifty,
              percent_aboveseventy, percent_aboveninety]
    toptitle = "Percent identity in markers"
    println(toptitle, ":\n", UnicodePlots.barplot(labels, values))
    Plots.pie(labels, values, palette = :Greens_8, title = toptitle, linecolor = "black")
    Plots.savefig(fileout)
    return fileout
end

