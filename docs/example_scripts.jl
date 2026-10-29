module ExampleScripts

using Markdown

function generate(source_dir, output_dir)
        mkpath(output_dir)
        scripts = String[]
        for source in sort(readdir(source_dir; join=true))
                endswith(source, ".md") || continue
                blocks = String[]
                contexts = Set{String}()
                for block in Markdown.parse(read(source, String)).content
                        block isa Markdown.Code || continue
                        header = split(block.language)
                        isempty(header) && continue
                        kind = first(header)
                        if kind in ("@example", "@setup")
                                context = strip(first(split(join(header[2:end], " "), ';')))
                                isempty(context) && error("Unnamed example context in $source")
                                push!(contexts, context)
                                push!(blocks, block.code)
                        elseif kind == "julia"
                                push!(blocks, block.code)
                        elseif kind in ("@repl", "jldoctest")
                                error("Unsupported executable block $kind in $source")
                        end
                end
                length(contexts) > 1 && error("Multiple example contexts in $source")
                isempty(blocks) && error("No Julia examples in $source")
                script = joinpath(output_dir, splitext(basename(source))[1] * ".jl")
                notice = "# This script is automatically generated without explanatory comments.\n" *
                         "# The documentation for this script exists at https://rightleftspin.github.io/LiNCEGE.jl/dev/examples/$(splitext(basename(source))[1]).html\n\n"
                write(script, notice * "output_dir = joinpath(@__DIR__, \"outputs\")\nmkpath(output_dir)\n\n" * join(blocks, "\n\n") * "\n")
                push!(scripts, script)
        end
        return scripts
end

function run_examples(scripts, project_dir)
        for script in scripts
                @info "Running example" script
                mktempdir() do directory
                        command = `$(Base.julia_cmd()) --startup-file=no --project=$project_dir $script`
                        run(addenv(Cmd(command; dir=directory), "GKSwstype" => "100"))
                end
        end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
        ExampleScripts.generate(joinpath(@__DIR__, "src", "examples"), joinpath(@__DIR__, "..", "examples"))
end
