from glycresoft_ms2_classification.apps.theoretical_search_space import search_space_builder


def main(ms1_results_file, digest_file, site_list_file):
    digest = search_space_builder.parse_digest(digest_file)
    task = search_space_builder.ExactSearchSpace(
        ms1_results_file, enzyme_info=digest.enzyme,
        constant_modifications=digest.constant_modifications,
        variable_modifications=digest.variable_modifications,
        site_list=site_list_file)
    task.run()


if __name__ == '__main__':
    import sys
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("ms1_results_file")
    app.add_argument("digest_file")
    app.add_argument("site_list_file")
    args = app.parse_args()
    main(args.ms1_results_file, args.digest_file, args.site_list_file)
