using DelimitedFiles


open("C:\\Users\\warner\\Downloads\\Lib80x\\xs_dir.txt") do f
    line = readline(f)
    print(typeof(f))
    print(line)
end
