    data["sw_diff_group1"] = data["sw_diff_group1"].strip(",").split(",")
    data["sw_diff_group2"] = data["sw_diff_group2"].strip(",").split(",")


   if data["plottype"] == "sw_diff":
        curr_time = time.time()
        filename = "SW_Diff_{}".format(curr_time)
        outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename),
                       "w")
        outfile.write(
            "Gene,Tran,Window_start,Grp1_count,Grp2_count,Ratio,Difference\n")
        transcripts = transcripts_org[[  # TODO: I might need to create copy
            "transcript", "cds_start", "cds_stop", "length,gene"
        ]]
        trandict = table2dict(transcripts, "transcript")
        gene_dict = table2dict(transcripts, "transcript")

        def group_profile_dict(sw_diff_group, grp_file_paths_dict):
            grp_profile_dict = {}
            for file_id in sw_diff_group:
                # printgrp1_file_paths_dict["riboseq"]
                sqlite_db = SqliteDict(
                    f"{grp_file_paths_dict['riboseq'][int(file_id)]}",
                    autocommit=False,
                    decode=my_decoder,
                )
                offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                # print"offsets", offsets
                for transcript in trandict.keys():
                    if transcript not in grp_profile_dict:
                        grp_profile_dict[transcript] = {}
                    try:
                        counts = sqlite_db[transcript]
                    except Exception:
                        continue
                    subprofile = build_profile(counts, offsets, False)
                    for pos in subprofile:
                        try:
                            grp_profile_dict[transcript][pos] += subprofile[
                                pos]
                        except Exception:
                            grp_profile_dict[transcript][pos] = subprofile[pos]
                sqlite_db.close()

        grp1_file_paths_dict = fetch_file_paths(data["sw_diff_group1"],
                                                organism)
        # printgrp1_file_paths_dict
        grp2_file_paths_dict = fetch_file_paths(data["sw_diff_group2"],
                                                organism)
        grp1_profile_dict = group_profile_dict(data["sw_diff_group1"],
                                               grp1_file_paths_dict)
        grp2_profile_dict = group_profile_dict(data["sw_diff_group2"],
                                               grp2_file_paths_dict)

        for tran in grp1_profile_dict:
            if tran in grp2_profile_dict:
                tranlen = trandict[tran]["length"]
                startpoint = 0
                endpoint = tranlen
                if data["custom_search_region"] != "whole_gene":
                    try:
                        cds_start = int(trandict[tran]["cds_start"])
                        cds_stop = int(trandict[tran]["cds_stop"])
                    except Exception:
                        continue
                    if data["custom_search_region"] == "five_leader":
                        startpoint = 0
                        endpoint = cds_start
                    elif data["custom_search_region"] == "cds":
                        startpoint = cds_start
                        endpoint = cds_stop
                    elif data["custom_search_region"] == "three_trailer":
                        startpoint = cds_stop
                        endpoint = tranlen
                for i in range(startpoint, endpoint,
                               data["sw_diff_step_size"]):
                    grp1_count = 1.0
                    grp2_count = 1.0
                    for x in range(i, i + data["sw_diff_window_size"]):
                        if x in grp1_profile_dict[tran]:
                            grp1_count += grp1_profile_dict[tran][x]
                        if x in grp2_profile_dict[tran]:
                            grp2_count += grp2_profile_dict[tran][x]
                    diff = abs(grp1_count - grp2_count)
                    if diff > data["sw_diff_min_diffsa"]:
                        ratio = grp1_count / grp2_count
                        gene = gene_dict[tran]
                        outfile.write(
                            f"{gene},{tran},{i},{grp1_count},{grp2_count},"
                            f"{ratio},{grp1_count - grp2_count}\n")

        return filename

    if data["plottype"] == "bulk_dl":
        curr_time = time.time()
        traninfo_dict = {}
        if data["region"] == "custom":
            transcrips = transcrips_org.loc[transcrips_org.transcript.isin(
                te_tranlist)]
        elif data["region"] == "principal":
            transcrips = transcrips_org.loc[transcrips_org.principle]
        result = transcripts[[
            "transcript", "gene", "cds_start", "cds_stop", "length"
        ]]
        for row in result:
            traninfo_dict[row[0]] = {
                "gene": row[1],
                "cds_start": row[2],
                "cds_stop": row[3],
                "length": int(row[4]),
            }
        identifiers = []
        filepaths = []
        map_types = ["unambig"]
        if data["readlen_ambig"]:
            map_types.append("ambig")
        for transcript in traninfo_dict:
            outfile = open(
                "{}/static/tmp/{}_{}".format(config.SCRIPT_LOC, transcript,
                                             curr_time),
                "w",
            )
            filepaths.append("{}_{}".format(transcript, curr_time))
            count_dict = {}
            for filetype in file_paths_dict:
                for file_id in file_paths_dict[filetype]:
                    filepath = file_paths_dict[filetype][file_id]
                    files = get_table("files")
                    files = files.loc[
                        files.file_id == file_id,
                        ["files_name", "file_description"]].iloc[0]
                    identifier = "{}_({})".format(files.file_name,
                                                  files.file_description)
                    identifiers.append(identifier)
                    count_dict[identifier] = {}
                    for i in range(0, traninfo_dict[transcript]["length"]):
                        count_dict[identifier][i] = 0
                    if os.path.isfile(filepath):
                        sqlite_db = SqliteDict(f"{filepath}",
                                               autocommit=False,
                                               decode=my_decoder)
                        offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                    else:
                        return (
                            f"File not found: {filepath}, please report this to"
                            " tripsvizsite@gmail.com or via the contact page.")
                    if transcript not in sqlite_db:
                        continue
                    for map_type in map_types:
                        counts = sqlite_db[transcript][map_type]
                        for readlen in counts:
                            if readlen > data["minreadlen"] and readlen < data[
                                    "maxreadlen"]:
                                if readlen in offsets:
                                    offset = offsets[readlen]
                                else:
                                    offset = 15
                                for pos in counts[readlen]:
                                    if apply_offset:
                                        try:
                                            count_dict[identifier][
                                                pos +
                                                offset] += counts[readlen][pos]
                                        except Exception:
                                            pass
                                    else:
                                        count_dict[identifier][pos] += counts[
                                            readlen][pos]
            if traninfo_dict[transcript]["cds_start"]:
                outfile.write("#{}_{}_{}_{}\n".format(
                    transcript,
                    traninfo_dict[transcript]["gene"],
                    traninfo_dict[transcript]["cds_start"],
                    traninfo_dict[transcript]["cds_stop"],
                ))
            else:
                outfile.write("#{}_{}_noncoding\n".format(
                    transcript, traninfo_dict[transcript]["gene"]))
            outfile.write("Pos,")
            for identifier in identifiers:
                outfile.write("{},".format(identifier))
            outfile.write("Total\n")
            for i in range(0, traninfo_dict[transcript]["length"]):
                outfile.write(str(i) + ",")
                pos_total = 0
                for identifier in identifiers:
                    outfile.write("{},".format(count_dict[identifier][i]))
                    pos_total += count_dict[identifier][i]
                outfile.write("{}\n".format(pos_total))

        subprocess.call(
            "tar -C {1}/static/tmp/ -czvf {1}/static/tmp/bulk_dl_{0}.tar.gz {2}"
            .format(
                curr_time,
                config.SCRIPT_LOC,
                str(filepaths).strip("[]").replace(",", "").replace("'", ""),
            ),
            shell=True,
        )
        return ("<div style='padding-left: 55px;padding-top: 22px;'><a href="
                f"'https://trips.ucc.ie/static/tmp/bulk_dl_{curr_time}.tar.gz'"
                " target='_blank' >"
                "<button class='button centerbutton' type='submit'>"
                "<b>Download result</b></button></a> </div>")

        if plottype == "mismatch_pos":
            master_dict = {}
            for filetype in file_paths_dict:
                for file_id in file_paths_dict[filetype]:
                    filepath = file_paths_dict[filetype][file_id]
                    if os.path.isfile(filepath):
                        sqlite_db = SqliteDict(filepath, autocommit=False)
                    else:
                        connection.close()
                        return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath))

                    if "global_mismatches" not in sqlite_db:
                        connection.close()
                        return ("No mismatch data for this file, please report this to tripsvizsite@gmail.com or via the contact page.")
                    else:
                        mismatches = sqlite_db["global_mismatches"]
                    sqlite_db.close()
                    for readlen in mismatches:
                        if readlen < mismatch_minreadlen or readlen > mismatch_maxreadlen:
                            continue
                        for pos in mismatches[readlen]:
                            if pos in master_dict:
                                master_dict[pos] += mismatches[readlen][pos]
                            else:
                                master_dict[pos] = mismatches[readlen][pos]
            title = "Mismatch positions"
            connection.close()

            return metainfo_plots.mismatch_pos(master_dict, title, short_code, background_col, readlength_col, title_size, axis_label_size, subheading_size, marker_size)

        elif plottype == "single_tran_de":
            if single_tran_de_range1.lower() in ["cds", "trailer", "leader"] or single_tran_de_range2.lower() in ["cds", "trailer", "leader"]:
                traninfo_connection = sqlite3.connect(
                    "/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{1}.sqlite".format(organism, transcriptome))
                traninfo_cursor = traninfo_connection.cursor()
                traninfo_cursor.execute("SELECT transcript,cds_start,cds_stop,length FROM transcripts WHERE transcript = '{}';".format(
                    single_tran_de_transcript))
                result = traninfo_cursor.fetchall()
                cds_start = int(result[0][1])
                cds_stop = int(result[0][2])
                length = int(result[0][3])
                if single_tran_de_range1.lower() == "cds":
                    range1 = [cds_start, cds_stop]
                if single_tran_de_range1.lower() == "leader":
                    range1 = [0, cds_start]
                if single_tran_de_range1.lower() == "trailer":
                    range1 = [cds_stop, length]
                if single_tran_de_range2.lower() == "cds":
                    range2 = [cds_start, cds_stop]
                if single_tran_de_range2.lower() == "leader":
                    range2 = [0, cds_start]
                if single_tran_de_range2.lower() == "trailer":
                    range2 = [cds_stop, length]
            else:
                range1 = single_tran_de_range1.split("_")
                # range1_kbp = (float(int(range1[1]) - int(range1[0])))/1000
                range2 = single_tran_de_range2.split("_")
            range1_len = int(range1[1])-int(range1[0])
            range2_len = int(range2[1])-int(range2[0])
            filename = "Single_tran_ratio_{}_{}_{}_{}_{}".format(
                single_tran_de_transcript, range1[0], range1[1], range2[0], range2[1])
            outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC, filename), "w")
            outfile.write(
                "Tran,Study,Range1_count, Range2_count, Norm_range1_count, Norm_range2_count,((Norm_range2_count/Norm_range1_count)*100)\n")

            master_list = []
            master_dict = {}
            for filetype in file_paths_dict:
                for file_id in file_paths_dict[filetype]:
                    range1_count = 0
                    range2_count = 0

                    filepath = file_paths_dict[filetype][file_id]
                    filename = filepath.split("/")[-1]
                    cursor.execute(
                        "SELECT file_description from files WHERE file_id = {};".format(file_id))
                    result = (cursor.fetchone())
                    file_desc = result[0]
                    study = filepath.split("/")[-2]
                    if study not in master_dict:
                        master_dict[study] = {"range1_count": 1.0, "range2_count": 1.0}
                    if os.path.isfile(filepath):
                        # Add the counts to the profile
                        sqlite_db = SqliteDict(filepath, autocommit=False)
                        if "offsets" in sqlite_db:
                            offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                        else:
                            offsets = {}
                        profile = {}
                        if single_tran_de_transcript in sqlite_db:
                            sqlite_db_tran = sqlite_db[single_tran_de_transcript]
                            for readlen in sqlite_db_tran["unambig"]:
                                if readlen in offsets:
                                    offset = offsets[readlen]
                                else:
                                    offset = 15
                                for pos in sqlite_db_tran["unambig"][readlen]:
                                    count = sqlite_db_tran["unambig"][readlen][pos]
                                    offset_pos = offset+pos
                                    if offset_pos not in profile:
                                        profile[offset_pos] = 0
                                    profile[offset_pos] += count
                        for x in range(int(range1[0]), int(range1[1])):
                            if x in profile:
                                range1_count += profile[x]
                        for x in range(int(range2[0]), int(range2[1])):
                            if x in profile:
                                range2_count += profile[x]
                    master_dict[study]["range1_count"] += range1_count
                    master_dict[study]["range2_count"] += range2_count
                    master_list.append((file_id, filename, range1_count+1.0,
                                       range2_count+1.0, file_desc, study))
            study_master_list = []
            for study in master_dict:
                range1_count = master_dict[study]["range1_count"]
                range2_count = master_dict[study]["range2_count"]
                study_master_list.append((0, study, range1_count, range2_count))
                norm_range1_count = range1_count/range1_len
                norm_range2_count = range2_count/range2_len
                per = (norm_range2_count/(norm_range1_count+1))*100
                outfile.write("{},{},{},{},{},{},{}\n".format(single_tran_de_transcript, study,
                              range1_count, range2_count, norm_range1_count, norm_range2_count, per))
            sorted_master_list = sorted(master_list, key=lambda x: x[1])
            connection.close()

            return metainfo_plots.single_tran_de(single_tran_de_transcript, sorted_master_list, study_master_list, organism, transcriptome, single_tran_de_study)

        elif plottype == "codon_usage":
            traninfo_connection = sqlite3.connect(
                "/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{1}.sqlite".format(organism, transcriptome))
            traninfo_cursor = traninfo_connection.cursor()
            codon_dict = {}
            principal_transcripts = {}
            traninfo_cursor.execute(
                "SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE principal = 1;")
            result = traninfo_cursor.fetchall()
            for row in result:
                if row[2] != "None" and row[2] != "" and row[2] != None:
                    principal_transcripts[str(row[0])] = {"seq": str(
                        row[1]), "cds_start": int(row[2]), "cds_stop": int(row[3])}

            if file_paths_dict["riboseq"] == {} and file_paths_dict["rnaseq"] == {}:
                flash("Error no files selected")
                return ("Error no files selected")
            all_values = []
            offset_dict = {}
            for file_id in file_paths_dict["riboseq"]:
                sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
                try:
                    offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                    offset_dict[file_id] = offsets
                except:
                    offset_dict[file_id] = {}
                sqlite_db.close()
            tran_count = 0

            for file_id in file_paths_dict["riboseq"]:
                sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
                if "codon_usage_dict2" in sqlite_db:
                    codon_usage_dict = sqlite_db["codon_usage_dict"]
                    for codon in codon_usage_dict:
                        if codon not in codon_dict:
                            codon_dict[codon] = {"ribo_count": 0, "codon_count": 0.0}
                        codon_dict[codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
                        codon_dict[codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
                else:
                    # codon_dict is the main dict that holds counts from all files, codon_usage_dict is file specific and will be saved for quick access later.
                    codon_usage_dict = {}
                    offsets = offset_dict[file_id]
                    # print"old offsets", offsets
                    poffsets = {}
                    for offset in offsets:
                        value = offsets[offset]
                        new_value = value - 3
                        poffsets[offset] = new_value
                    # print"new offsets", poffsets
                    for transcript in principal_transcripts:
                        tran_count += 1
                        profile = {}
                        if transcript not in sqlite_db:
                            continue

                        subprofile = build_profile(sqlite_db[transcript], poffsets, "unambig")
                        for pos in subprofile:
                            if pos not in profile:
                                profile[pos] = 0
                            profile[pos] += subprofile[pos]
                        seq = principal_transcripts[transcript]["seq"]
                        for i in range(principal_transcripts[transcript]["cds_start"]-1, principal_transcripts[transcript]["cds_start"]+30):
                            codon = seq[i:i+3]
                            if len(codon) != 3:
                                continue
                            count = 0
                            if i in profile:
                                count += profile[i]
                            if codon not in codon_dict:
                                codon_dict[codon] = {"ribo_count": 0, "codon_count": 0.0}
                            if codon not in codon_usage_dict:
                                codon_usage_dict[codon] = {"ribo_count": 0, "codon_count": 0.0}
                            codon_usage_dict[codon]["ribo_count"] += count
                            codon_usage_dict[codon]["codon_count"] += 1
                    for codon in codon_usage_dict:
                        codon_dict[codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
                        codon_dict[codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
                    sqlite_db["codon_usage_dict"] = codon_usage_dict
                    sqlite_db.commit()
                sqlite_db.close()
            connection.close()
            return metainfo_plots.codon_usage(codon_dict, short_code, str(title_size)+"pt", str(axis_label_size)+"pt", str(marker_size)+"pt")

        elif plottype == "diff_codon_usage":
            traninfo_connection = sqlite3.connect(
                "/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{1}.sqlite".format(organism, transcriptome))
            traninfo_cursor = traninfo_connection.cursor()
            codon_dict_cond = {"condition1": {}, "condition2": {}}
            diff_codon_dict = {}
            principal_transcripts = {}
            condition_totals = {"condition1": 0, "condition2": 0}
            traninfo_cursor.execute(
                "SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE principal = 1;")
            result = traninfo_cursor.fetchall()
            for row in result:
                if row[2] != "None" and row[2] != "" and row[2] != None:
                    principal_transcripts[str(row[0])] = {"seq": str(
                        row[1]), "cds_start": int(row[2]), "cds_stop": int(row[3])}

            if file_paths_dict["riboseq"] == {} and file_paths_dict["rnaseq"] == {}:
                flash("Error no files selected")
                connection.close()
                return ("Error no files selected")

            condition_dict = {"condition1": [file_paths_dict["riboseq"].keys(
            )[0]], "condition2": [file_paths_dict["riboseq"].keys()[1]]}

            all_values = []
            offset_dict = {}
            for condition in condition_dict:
                for file_id in condition_dict[condition]:
                    sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
                    try:
                        offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
                        offset_dict[file_id] = offsets
                    except:
                        offset_dict[file_id] = {}
                    sqlite_db.close()
                tran_count = 0

                for file_id in condition_dict[condition]:
                    sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
                    if "codon_usage_dict" in sqlite_db:
                        codon_usage_dict = sqlite_db["codon_usage_dict"]
                        for codon in codon_usage_dict:
                            if codon not in codon_dict_cond[condition]:
                                codon_dict_cond[condition][codon] = {
                                    "ribo_count": 0, "codon_count": 0.0}
                            codon_dict_cond[condition][codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
                            condition_totals[condition] += codon_usage_dict[codon]["ribo_count"]
                            codon_dict_cond[condition][codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
                    else:
                        # codon_dict_cond is the main dict that holds counts from all files, codon_usage_dict is file specific and will be saved for quick access later.
                        codon_usage_dict = {}
                        for transcript in principal_transcripts:
                            tran_count += 1
                            profile = {}
                            if transcript not in sqlite_db:
                                continue
                            offsets = offset_dict[file_id]
                            subprofile = build_profile(sqlite_db[transcript], offsets, "unambig")
                            for pos in subprofile:
                                if pos not in profile:
                                    profile[pos] = 0
                                profile[pos] += subprofile[pos]
                            seq = principal_transcripts[transcript]["seq"]
                            for i in range(principal_transcripts[transcript]["cds_start"], principal_transcripts[transcript]["cds_stop"], 3):
                                codon = seq[i:i+3]
                                if len(codon) != 3:
                                    continue
                                count = 0
                                if i in profile:
                                    count += profile[i]
                                if i+1 in profile:
                                    count += profile[i+1]
                                if i+2 in profile:
                                    count += profile[i+2]
                                if codon not in codon_dict_cond[condition]:
                                    codon_dict_cond[condition][codon] = {
                                        "ribo_count": 0, "codon_count": 0.0}
                                # codon_dict_cond[codon]["ribo_count"] += count
                                # codon_dict_cond[codon]["codon_count"] += 1

                                if codon not in codon_usage_dict:
                                    codon_usage_dict[codon] = {"ribo_count": 0, "codon_count": 0.0}
                                codon_usage_dict[codon]["ribo_count"] += count
                                codon_usage_dict[codon]["codon_count"] += 1
                        for codon in codon_usage_dict:
                            codon_dict_cond[condition][codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
                            condition_totals[condition] += codon_usage_dict[codon]["ribo_count"]
                            codon_dict_cond[condition][codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
                        sqlite_db["codon_usage_dict"] = codon_usage_dict
                        sqlite_db.commit()
                    sqlite_db.close()
            factor_diff = float(
                condition_totals["condition1"])/float(condition_totals["condition2"])
            for codon in codon_dict_cond["condition1"]:
                count1 = codon_dict_cond["condition1"][codon]["ribo_count"]
                count2 = codon_dict_cond["condition2"][codon]["ribo_count"]*factor_diff
                diff = count1-count2
                diff_codon_dict[codon] = {
                    "ribo_count": diff, "codon_count": codon_dict_cond["condition1"][codon]["codon_count"]}
            connection.close()
            return metainfo_plots.codon_usage(diff_codon_dict, short_code, str(title_size)+"pt", str(axis_label_size)+"pt", str(marker_size)+"pt")
	elif plottype == "tran_corr":
		master_list = []
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				tran1_count = 1.001
				tran2_count = 1.001
				filepath = file_paths_dict[filetype][file_id]
				filename = filepath.split("/")[-1]
				study = filepath.split("/")[-2]
				if filename not in master_dict:
					master_dict[filename] = {"tran1_count": 0,
											"tran2_count": 0}
				if os.path.isfile(filepath):
					# Add the counts to the profile
					sqlite_db = SqliteDict(filepath, autocommit=False)
					if "offsets" in sqlite_db:
						offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
					else:
						offsets = {}
					profile = {}
					# TRAN1
					if tran_corr_transcript1 in sqlite_db:
						sqlite_db_tran = sqlite_db[tran_corr_transcript1]
						for readlen in sqlite_db_tran["unambig"]:
							if readlen in offsets:
								offset = offsets[readlen]
							else:
								offset = 15
							for pos in sqlite_db_tran["unambig"][readlen]:
								count = sqlite_db_tran["unambig"][readlen][pos]
								offset_pos = offset+pos
								if offset_pos not in profile:
									profile[offset_pos] = 0
								profile[offset_pos] += count
					for pos in profile:
						tran1_count += profile[pos]
					# TRAN2
					profile = {}
					if tran_corr_transcript2 in sqlite_db:
						sqlite_db_tran = sqlite_db[tran_corr_transcript2]
						for readlen in sqlite_db_tran["unambig"]:
							if readlen in offsets:
								offset = offsets[readlen]
							else:
								offset = 15
							for pos in sqlite_db_tran["unambig"][readlen]:
								count = sqlite_db_tran["unambig"][readlen][pos]
								offset_pos = offset+pos
								if offset_pos not in profile:
									profile[offset_pos] = 0
								profile[offset_pos] += count
					for pos in profile:
						tran2_count += profile[pos]
				master_dict[filename]["tran1_count"] += tran1_count
				master_dict[filename]["tran2_count"] += tran2_count
				master_list.append((file_id, filename, log(
				    tran1_count, 2), log(tran2_count, 2), study))
		sorted_master_list = sorted(master_list, key=lambda x: x[2])
		connection.close()
		return metainfo_plots.tran_corr(tran_corr_transcript1, tran_corr_transcript2, sorted_master_list, organism, transcriptome)

	elif plottype == "mismatches":
		positive_hits = 0
		result_list = []
		file_string = ""
		total_trans = 0
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				file_string += "{},".format(file_id)
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(
		    organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			traninfo_dict = SqliteDict("{0}/{1}/{2}/{2}.sqlite".format(
			    config.SCRIPT_LOC, config.ANNOTATION_DIR, organism), autocommit=False)
		else:
			traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(
			    config.UPLOADS_DIR, owner, organism, transcriptome), autocommit=False)
		if organism == "homo_sapiens" or organism == "homo_sapiens_polio":
			longest_tran_db = SqliteDict("{0}/{1}/homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite".format(
			    config.SCRIPT_LOC, config.ANNOTATION_DIR), autocommit=False)
			longest_tran_list = longest_tran_db["transcripts"]
			longest_tran_db.close()
		else:
			longest_tran_list = traninfo_dict.keys()

		if mismatch_agg == True:
			for transcript in longest_tran_list:
				total_trans += 1

				cds_start = traninfo_dict[transcript]["cds_start"]
				cds_stop = traninfo_dict[transcript]["cds_stop"]
				tranlen = traninfo_dict[transcript]["length"]
				if mismatch_region == "all":
					minpos = 0
					maxpos = tranlen
				elif mismatch_region == "fiveleader":
					minpos = 0
					maxpos = cds_start
				elif mismatch_region == "cds":
					minpos = cds_start
					maxpos = cds_stop
				elif mismatch_region == "threetrailer":
					minpos = cds_stop
					maxpos = tranlen
				elif mismatch_region == "cds_start":
					minpos = cds_start-0
					maxpos = cds_start+3
				elif mismatch_region == "cds_stop":
					minpos = cds_stop-1
					maxpos = cds_stop+2
				if positive_hits > mismatch_maxhit:
					break
				sequence = traninfo_dict[transcript]["seq"]
				profile = {}
				mismatch_profile = {"A": {}, "T": {}, "G": {}, "C": {}}
				for filetype in file_paths_dict:
					for file_id in file_paths_dict[filetype]:

						filepath = file_paths_dict[filetype][file_id]
						if os.path.isfile(filepath):

							# Add the counts to the profile
							sqlite_db = SqliteDict(filepath, autocommit=False)
							if transcript in sqlite_db:
								sqlite_db_tran = sqlite_db[transcript]
								for readlen in sqlite_db_tran["unambig"]:
									for pos in sqlite_db_tran["unambig"][readlen]:
										count = sqlite_db_tran["unambig"][readlen][pos]
										for x in range(pos, pos+readlen):
											if x not in profile:
												profile[x] = 0
											profile[x] += count

								# Add the mismatch to the profile
								sqlite_db_seqvar = dict(sqlite_db[transcript]["seq"])
								for pos in sqlite_db_seqvar:
									# convert to one based
									fixed_pos = pos+1
									for char in sqlite_db_seqvar[pos]:
										if char != "N":
											if fixed_pos not in mismatch_profile[char]:
												mismatch_profile[char][fixed_pos] = 0
											count = sqlite_db_seqvar[pos][char]
											mismatch_profile[char][fixed_pos] += count
						else:
							connection.close()
							return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath))
						sqlite_db.close()

				for pos in profile:
					if pos < minpos or pos > maxpos:
						continue
					if profile[pos] > mismatch_minreadcount:
						for char in mismatch_profile:
							if pos in mismatch_profile[char]:
								mismatch_count = mismatch_profile[char][pos]
								per = (float(mismatch_count)/float(profile[pos]))*100
								if per > 100:
									per = 100
								transition = "{}->{}".format(sequence[pos-1], char)
								if per > mismatch_minper and per < mismatch_maxper:
									# return "Mismatch in transcript {} at position {}, read count:{}, mismatch count {}, transition {}".format(transcript, pos, profile[pos], mismatch_count, transition)
									trips_link = '<a href="https://trips.ucc.ie/'+organism+'/'+transcriptome+'/interactive_plot/?&hili=' + \
									    str(pos-10)+'_'+str(pos+10)+'&tran='+transcript+'&cov=T&nuc=T&files=' + \
									        file_string+'" target="_blank_" >View on trips-viz</a>'
									result_list.append([transcript, pos, profile[pos], mismatch_count,
									                   transition, int(per), trips_link, "Aggregate"])
									positive_hits += 1
		else:
			for transcript in longest_tran_list:
				total_trans += 1
				cds_start = traninfo_dict[transcript]["cds_start"]
				cds_stop = traninfo_dict[transcript]["cds_stop"]
				tranlen = traninfo_dict[transcript]["length"]
				if mismatch_region == "all":
					minpos = 0
					maxpos = tranlen
				elif mismatch_region == "fiveleader":
					minpos = 0
					maxpos = cds_start
				elif mismatch_region == "cds":
					minpos = cds_start
					maxpos = cds_stop
				elif mismatch_region == "threetrailer":
					minpos = cds_stop
					maxpos = tranlen
				elif mismatch_region == "cds_start":
					minpos = cds_start-0
					maxpos = cds_start+2
				elif mismatch_region == "cds_stop":
					minpos = cds_stop-0
					maxpos = cds_stop+2

				if positive_hits > mismatch_maxhit:
					break
				sequence = traninfo_dict[transcript]["seq"]
				file_string = ""
				for filetype in file_paths_dict:
					for file_id in file_paths_dict[filetype]:
						file_string += "{},".format(file_id)
						filepath = file_paths_dict[filetype][file_id]
						if os.path.isfile(filepath):

							# Add the counts to the profile
							sqlite_db = SqliteDict(filepath, autocommit=False)
							profile = {}
							mismatch_profile = {"A": {}, "T": {}, "G": {}, "C": {}}
							if transcript in sqlite_db:
								sqlite_db_tran = sqlite_db[transcript]
								for readlen in sqlite_db_tran["unambig"]:
									for pos in sqlite_db_tran["unambig"][readlen]:
										count = sqlite_db_tran["unambig"][readlen][pos]
										for x in range(pos, pos+readlen):
											if x not in profile:
												profile[x] = 0
											profile[x] += count

								# Add the mismatch to the profile
								sqlite_db_seqvar = dict(sqlite_db[transcript]["seq"])
								for pos in sqlite_db_seqvar:
									# convert to one based
									fixed_pos = pos+1
									for char in sqlite_db_seqvar[pos]:
										if char != "N":
											if fixed_pos not in mismatch_profile[char]:
												mismatch_profile[char][fixed_pos] = 0
											count = sqlite_db_seqvar[pos][char]
											mismatch_profile[char][fixed_pos] += count
						else:
							connection.close()
							return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath))
						sqlite_db.close()

						for pos in profile:
							if pos < minpos or pos > maxpos:
								continue
							if profile[pos] > mismatch_minreadcount:
								for char in mismatch_profile:
									if pos in mismatch_profile[char]:
										mismatch_count = mismatch_profile[char][pos]
										per = (float(mismatch_count)/float(profile[pos]))*100
										if per > 100:
											per = 100
										transition = "{}->{}".format(sequence[pos-1], char)
										if per > mismatch_minper and per < mismatch_maxper:
											# return "Mismatch in transcript {} at position {}, read count:{}, mismatch count {}, transition {}".format(transcript, pos, profile[pos], mismatch_count, transition)
											trips_link = '<a href="https://trips.ucc.ie/'+organism+'/'+transcriptome+'/interactive_plot/?&hili=' + \
											    str(pos-10)+'_'+str(pos+10)+'&tran='+transcript+'&cov=T&nuc=T&files=' + \
											        file_string+'" target="_blank_" >View on trips-viz</a>'
											result_list.append([transcript, pos, profile[pos], mismatch_count, transition, int(
											    per), trips_link, (filepath.split("/"))[-1]])
											positive_hits += 1

		table_str = "<table class='prediction_table hover'><tr><th>Transcript</th><th>Position</th><th>Filename</th><th>Read Count</th><th>Mismatch Count</th><th>Percentage</th><th>Transition</th><th>View</th></td>"
		for item in result_list:
			table_str += "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
			    item[0], item[1], item[7], item[2], item[3], item[5], item[4], item[6])
		table_str += "</table>"
		connection.close()
		return table_str




	elif plottype == "explore_offsets":
		readlen_dict = {}
		traninfo_dict = SqliteDict("{0}/{1}/{2}/{2}.{3}.sqlite".format(
		    config.SCRIPT_LOC, config.ANNOTATION_DIR, organism, transcriptome), autocommit=False)
		tranlist = traninfo_dict.keys()[:10000]
		labels = []
		f0_counts = []
		f1_counts = []
		f2_counts = []

		# For the first file in selected files
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict = dict(sqlite_db)
					sqlite_db.close()
				else:
					connection.close()
					return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page".format(filepath))
				tranlist = opendict.keys()

				# For each readlength we display on the final graph
				for readlen in range(25, 36):
					try:
						chosen_offset = opendict["offsets"]["fiveprime"]["offsets"][readlen]
					except:
						chosen_offset = 15
					labels.append("")
					labels.append("{}_{}".format(readlen, chosen_offset))
					labels.append("")
					for offset in [chosen_offset-1, chosen_offset, chosen_offset+1]:
						trancount = 0
						inframe_counts = 0
						outframe_counts = 0
						if readlen not in readlen_dict:
							readlen_dict[readlen] = {chosen_offset-1: [0, 0, 0],
							    chosen_offset: [0, 0, 0], chosen_offset+1: [0, 0, 0]}
						# for each transcript get the frame counts breakdown from the cds given this particular offset
						for tran in tranlist:
							if tran not in traninfo_dict:
								continue
							tempdict = dict(opendict[tran])
							trancount += 1
							if trancount > 5000:
								break

							if "cds_start" not in traninfo_dict[tran]:
								continue
							cds_start = traninfo_dict[tran]["cds_start"]
							cds_stop = traninfo_dict[tran]["cds_stop"]

							if cds_start == "NULL" or cds_stop == "NULL":
								continue
							if cds_start <= 1 or cds_stop <= 1:
								continue
							# to account for 0-based counts ,without this line the frame will be wrong
							cds_start += 1
							cds_frame = cds_start % 3
							# first walk through this entry in the opendict for only the readlength in question applying the relevant offset
							count_dict = {}

							if readlen in tempdict["unambig"]:
								for fiveprime_pos in tempdict["unambig"][readlen]:
									count = tempdict["unambig"][readlen][fiveprime_pos]
									new_pos = fiveprime_pos + offset
									count_dict[new_pos] = count
								for i in range(cds_start, cds_stop):
									frame = i % 3
									if i in count_dict:
										if frame == cds_frame:
											inframe_counts += count_dict[i]
										else:
											outframe_counts += count_dict[i]

						f0_counts.append(inframe_counts)
						f1_counts.append(outframe_counts)
						f2_counts.append(0)
		connection.close()
		return metainfo_plots.explore_offsets(f0_counts, f1_counts, f2_counts, labels, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size)


	elif plottype == "rust_dwell":
		traninfo_dict = SqliteDict("{0}/{1}/{2}/{2}.{3}.sqlite".format(
		    config.SCRIPT_LOC, config.ANNOTATION_DIR, organism, transcriptome), autocommit=False)
		longest_tran_db = SqliteDict(
		    "/home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite", autocommit=True)
		longest_tran_list = longest_tran_db["transcripts"]
		codon_count_dict = {"TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
		"TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
		"TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
		"TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
		"CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
		"CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
		"CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
		"CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
		"ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
		"ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
		"AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
		"AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
		"GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
		"GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
		"GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
		"GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict = dict(sqlite_db)
					sqlite_db.close()
				else:
					connection.close()
					return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath))
				# TODO CHANGE THIS SO THAT THE CODON COUNT DICT IS OUTSIDE THE FILEPATH FOR LOOP
				offsets = opendict["offsets"]["fiveprime"]["offsets"]
				position_dict = {}
				for transcript in longest_tran_list:
					if transcript in opendict:
						cds_start = traninfo_dict[transcript]["cds_start"]
						cds_stop = traninfo_dict[transcript]["cds_stop"]
						transeq = traninfo_dict[transcript]["seq"]
						for readlen in opendict[transcript]["unambig"]:
							offset = 15
							for pos in opendict[transcript]["unambig"][readlen]:
								a_site = (pos+offset)+1
								if a_site > cds_start+120 and a_site < cds_stop-60:
									codon = transeq[a_site:a_site+3]
									codon_count_dict[codon] += opendict[transcript]["unambig"][readlen][pos]
		connection.close()
		return metainfo_plots.rust_dwell(codon_count_dict, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size)

	elif plottype == "unmapped":
		return metainfo_plots.most_freq_unmapped(file_paths_dict, short_code)
	elif plottype == "contamination":
		count_dict = {}
		master_sequence = ""
		contaminant_file = open(
		    "/home/DATA/www/tripsviz/tripsviz/static/contaminants/mycoplasma.fa")
		contaminant_lines = contaminant_file.read()
		contaminant_split = contaminant_lines.split(">")
		for entry in contaminant_split[1:]:
			header = entry.split("\n")[0]
			sequence = "".join(entry.split("\n")[1:])
			master_sequence += sequence
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				filename = filepath.split("/")[-1]
				if filename not in count_dict:
					count_dict[filename] = {"count": 0, "coverage": [], "unique_reads": 0}
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					connection.close()
					return ("File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1]))
				if "frequent_unmapped_reads" not in sqlite_db:
					connection.close()
					return ("No unmapped reads data for {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1]))
				# unmapped reads list is a list of tuples of length 100, first item in tuple is a sequence second is a count
				unmapped_reads_list = sqlite_db["frequent_unmapped_reads"]
				sqlite_db.close()
				for tup in unmapped_reads_list:
					read = tup[0]
					count = tup[1]
					readlen = len(read)
					for x in range(0, len(master_sequence)-readlen):
						mismatches = 0
						for y in range(0, readlen):
							if master_sequence[x+y] != read[y]:
								mismatches += 1
							if mismatches > 2:
								break
						if mismatches <= 2:
							count_dict[filename]["unique_reads"] += 1
							count_dict[filename]["count"] += count
							for i in range(x, x+readlen):
								if i not in count_dict[filename]["coverage"]:
									count_dict[filename]["coverage"].append(i)
		master_seq_len = len(master_sequence)
		for filename in count_dict:
			coverage = float(
			    len(count_dict[filename]["coverage"]))/float(master_seq_len)
			coverage = round((coverage*100), 2)
			count_dict[filename]["coverage"] = coverage
		title = "Contamination counts ({})".format(short_code)
		top_reads = (sorted(count_dict.items(), key=operator.itemgetter(1)))
		html_table = "<h1><center>{}</center></h1>".format(title)
		html_table += """<table class="unmapped_table">
		<thead><tr><th>Filename</th><th>Counts</th><th>Unique reads</th><th>Percentage coverage</th></tr></thead>"""
		for tup in top_reads[::-1]:
			html_table += ("<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>".format(
			    tup[0], tup[1]["count"], tup[1]["unique_reads"], tup[1]["coverage"]))
		html_table += ("</table>")
		connection.close()
		return html_table
	elif plottype == "fastq_screen":
		html_filepath = ""
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if html_filepath == "":
					html_filepath = filepath.replace(".sqlite", "_lessrRNA_screen.html")
				else:
					connection.close()
					return ("Error: Only one dataset at a time can be selected for fastq screen")
		if os.path.isfile(html_filepath):
			openfile = open(html_filepath, "r")
			fastq_lines = openfile.readlines()
			# The base64 encoded png string in the header is too long for firefox, will work for one plot and then crash firefox, this is a hack to prevent that
			fixed_html = ""
			for line in fastq_lines:
				if "iVBORw0KGgoAAAANSUhEUgAAA4wAAAGVCAYAAAHC" in line:
					fixed_html += '<a style="float:left;" href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen" target="_blank"><img width="50%" height="50%" alt="FastQ Screen"src="/static/fastq_screen.png"</a>'
					continue
				else:
					fixed_html += (line)
			# Replace serves two functions here, first remove the padding line in the body tag as this affects the header bar and everything else on trips,
			# but removal means fastq screen logo is slightly off screen
			# Second remove the max-width line in the .container class, replace it with the padding line removed from the body tag, as this will now be specific to the container
			# and fix the fastq screen logo.
			fixed_html = str(fixed_html.replace("padding:0 20px 20px", "").replace("max-width:1200px;", "padding:0 20px 20px").replace("<html>", "").replace("</html>",
			                 "").replace("<body>", "").replace("</body>", "")).replace("<!DOCTYPE html>", "").replace("<head>", "").replace("</head>", "").replace("container", "container2")
			connection.close()
			return fixed_html
