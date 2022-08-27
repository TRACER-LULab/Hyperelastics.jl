bib_file = parse_file("../CITATIONS.bib")
get_model(s::String) = select(bib_file, "")
