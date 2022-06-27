module RunTests

println("Running tests...")
include("continuous_test.jl")
include("discrete_test.jl")
include("example_1_test.jl")
include("example_2_test.jl")
include("test_id.jl")
include("test_version.jl")
include("test_periodic.jl")

end # module
