using JuMP, GLPK

function readFile(filePath)
    s = open(filePath)

    vertices = parse(Int,split(readline(s))[2])
    paths = parse(Int,split(readline(s))[2])
    tmax = parse(Float32,split(readline(s))[2])

    prices = zeros(vertices)
    costs = zeros(vertices,vertices)

    xCoord = zeros(vertices)
    yCoord = zeros(vertices)
    for i in 1:vertices
       xCoord[i], yCoord[i], prices[i] = [parse(Float32,i) for i in split(readline(s))]
   end

    for i in 1:vertices-1
      for j in i+1:vertices
        costs[i,j] = ( (xCoord[i] - xCoord[j])^2 + (yCoord[i] - yCoord[j])^2 ) ^ 0.5
		costs[j,i] = costs[i,j]
      end
    end

    close(s)
    return (vertices,paths,tmax,prices,costs)
end

println("Starting")

dim, m, tmax, p, t = readFile(ARGS[1])

model = Model(with_optimizer(GLPK.Optimizer, tm_lim = 60000, msg_lev = GLPK.ON))

#variavel Yij: 1 se o arco (i,j) e utilizado, o caso contrario
@variable(model, y[i = 1:dim-1, j= 2:dim; i != j], Bin)
#variavel Fij: 0 <= Fij <= capacidade, fluxo que passa pelo arco (i,j)
@variable(model, 0 <= f[i = 1:dim-1, j= 2:dim; i != j] <= tmax)

#Funcao objetivo: SUM Pi * Yij
@expression(model, obj, sum(p[j]*y[i,j] for i in 1:dim-1 for j in 2:dim if i != j))
@objective(model, Max, obj)

#Restricao de numero de rotas: SUM Y1j <= M
@expression(model, originLH, sum(y[1,j] for j in 2:dim))
@constraint(model, origin, originLH <= m)

#Restricao de grau de um vertice: SUM Yij <= 1
@expression(model, grauINLH[i = 2:dim-1], sum(y[j,i] for j in 1:dim-1 if i !=j))
@constraint(model,grauIN[i = 2:dim-1], grauINLH[i] <= 1)

#Restricao de grau de um vertice: SUM Yij = SUM Yji
@expression(model, grauLH[i = 2:dim-1], sum(y[i,j] for j in 2:dim if i != j) - sum(y[j,i] for j in 1:dim-1 if i != j))
@constraint(model, grau[i = 2:dim-1], grauLH[i] == 0)

#Restricao de valor do fluxo de um arco: Fij <= capacidade * Yij
@constraint(model, flowValue[i = 1:dim-1, j = 2:dim; i != j], f[i,j] <= tmax*y[i,j])

#Restricao de sequencia de fluxo: SUM Fij - SUM Fji = SUM Yij * Tij
@expression(model, flowLH[i = 1:dim-1], sum(f[i,j] for j in 2:dim if i != j) - sum(f[j,i] for j in 1:dim-1 if i != j && i != 1 ))
@expression(model, flowRH[i = 1:dim-1], sum(y[i,j]*t[i,j] for j in 2:dim if i != j))
@constraint(model, flowSeq[i = 1:dim-1], flowLH[i] == flowRH[i])

open("modelo.lp", "w") do f
    print(f,model)
end

optimize!(model)

println("objective cost: ", objective_value(model))
arcs = [(i,j) for i in 1:dim-1 for j in 2:dim if i != j && value(y[i,j]) == 1]
sort!(arcs)
println(arcs)
flux = [(i,j,value(f[i,j])) for i in 1:dim-1 for j in 2:dim if i != j && value(f[i,j]) >= 0.1]
sort!(flux)
println(flux)

