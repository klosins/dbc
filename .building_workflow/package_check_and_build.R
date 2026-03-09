# Rebuild and load
devtools::load_all()

# Knit README
rmarkdown::render(
    "readme_code/README.Rmd",
    output_file = "../README.md",
    knit_root_dir = "."
)

# Check documentation, tests and spelling
devtools::document()
devtools::test()
devtools::check_man()
devtools::spell_check()

# Do local check (of build)
devtools::check()

# Check windows system
devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()

# Check macos
rhub::check(platform = "macos")

# Check list, ready for release?
usethis::use_release_issue()
