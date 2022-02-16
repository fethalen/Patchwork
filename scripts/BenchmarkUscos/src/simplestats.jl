using DataFrames
using Plots

function plot_querycover(
    results::DataFrame,
    outdir::AbstractString
)
    categories = [(0.0, 20.0), (21.0, 40.0), (41.0, 60.0), (61.0, 80.0), (81.0, 100.0)]
    labels = ["0.1-20%", "21-40%", "41-60%", "61-80%", "81-100%"]
    bincount = map(category ->
        count(x -> first(category) <= x <= last(category),
        results.query_cover),
        categories)
    toptitle = "Query coverage in markers"
    bar(labels, bincount, color = "green", title = toptitle, bar_width = 0.3,
        legend = false)
    figout = outdir * "/" * "query_coverage.png"
    savefig(figout)
    return figout
end

function plot_percentident(
    results::DataFrame,
    subjectcount::Int,
    outdir::AbstractString
)
    total = length(eachrow(results))
    @assert subjectcount > 0
    @assert subjectcount >= total
    aboveninety = length(eachrow(filter(df -> df.percent_identity > 90.0, results)))
    aboveseventy = length(eachrow(filter(df -> df.percent_identity > 70.0, results)))
    abovefifty = length(eachrow(filter(df -> df.percent_identity > 50.0, results)))
    abovethirty = length(eachrow(filter(df -> df.percent_identity > 30.0, results)))
    belowthirty = total - abovefifty
    labels = ["missing", "≤ 30%", "> 30%", ">50%", "> 70%", "> 90%"]

    missing = round(((subjectcount - total) / subjectcount) * 100.0, digits = 2)
    percent_aboveninety = round((aboveninety / subjectcount) * 100.0, digits = 2)
    percent_aboveseventy = round((aboveseventy / subjectcount) * 100.0, digits = 2)
    percent_abovefifty = round((abovefifty / subjectcount) * 100.0, digits = 2)
    percent_abovethirty = round((abovethirty / subjectcount) * 100.0, digits = 2)
    percent_belowthirty = round((belowthirty / subjectcount) * 100.0, digits = 2)
    values = [missing, percent_belowthirty, percent_abovethirty, percent_abovefifty,
              percent_aboveseventy, percent_aboveninety]
    # TODO: display percentage inside diagram
    toptitle = "Percent identity in markers"
    pie(labels, values, palette = :Greens_8, title = toptitle, linecolor = "black")
    figout = outdir * "/" * "percent_identity.png"
    savefig(figout)
    return figout
end
