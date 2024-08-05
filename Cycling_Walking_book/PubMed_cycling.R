### Exproer PubMed with R
### package: rentrez
### 2024/07/16 kizen sasaki

library(rentrez)
library(dplyr)
library(stringr)

# 検索クエリの設定
query <- "(cycling OR bicycling OR walking) AND (health OR fitness) 
  AND (commuting OR transportation) AND adult* 
  AND (\"2010\"[Date - Publication] : \"3000\"[Date - Publication])"

# Web Historyを使用して検索を実行
search_results <- entrez_search(db = "pubmed", 
                                term = query, use_history = TRUE)

# 総件数を取得
total_count <- search_results$count
cat("検索結果の総数:", total_count, "\n")

# バッチサイズを設定（PubMedの推奨は最大200）
batch_size <- 200

# 検索結果の詳細を取得（バッチ処理）
all_summaries <- list()
for(i in seq(1, total_count, batch_size)) {
  end <- min(i + batch_size - 1, total_count)
  tryCatch({
    batch_summaries <- entrez_summary(db = "pubmed", 
                                      web_history = search_results$web_history, 
                                      retstart = i - 1, 
                                      retmax = batch_size)
    all_summaries <- c(all_summaries, batch_summaries)
  }, error = function(e) {
    cat("Error occurred at batch", i, "to", end, ":", conditionMessage(e), "\n")
  })
  Sys.sleep(1)  # APIへの負荷を減らすために待機時間を1秒に設定
}

# NULLの要素を除去
all_summaries <- all_summaries[!sapply(all_summaries, is.null)]

# データフレームの作成
df <- do.call(rbind, lapply(all_summaries, function(x) {
  data.frame(
    PMID = ifelse(is.null(x$uid), NA, x$uid),
    Author = ifelse(is.null(x$sortfirstauthor), NA, x$sortfirstauthor),
    Year = ifelse(is.null(x$pubdate), NA, substr(x$pubdate, 1, 4)),
    Title = ifelse(is.null(x$title), NA, x$title),
    Abstract = ifelse(is.null(x$abstract), NA, x$abstract),
    stringsAsFactors = FALSE
  )
}))

# Study列の作成（Author Year の形式）
df$Study <- paste(df$Author, df$Year)

# CSVファイルとして保存
write.csv(df, "dfMeta_initial.csv", row.names = FALSE)

# 保存したファイルの行数を表示
cat("Saved", nrow(df), "records to dfMeta_initial.csv\n")
###### End of the code for all result ######

### CONSORT声明に従って、1回目のスクリーニングを行い、
### 除外基準に基づいて論文を選別する。除外する論文を特定し、
### その情報を"dfMeat.csv"ファイルに保存する。
### Excl1列に著者名と年 Oja, 2011 の形式で記入する。
# 必要なライブラリの読み込み
library(dplyr)
library(stringr)

df <- read.csv("dfMeta_cycling.csv", stringsAsFactors = FALSE)

# 除外基準の関数を定義
is_excluded <- function(title, abstract_col) {
  exclusion_keywords <- c("children", "adolescent", "review", "meta-analysis", 
                          "protocol", "case report", "letter", "comment")
  
  lower_title <- tolower(title)
  lower_abstract <- tolower(abstract_col)
  
  # タイトルに除外キーワードが含まれている場合のみ除外
  title_excluded <- any(sapply(exclusion_keywords, function(keyword) str_detect(lower_title, keyword)))
  
  # 抄録に特定のキーワードが含まれている場合は除外しない
  abstract_included <- grepl("(cycling|bicycling|walking)|(health|fitness|\"health-related quality of life\"|\"cardiorespiratory fitness\"|\"muscle strength\"|flexibility|balance)", lower_abstract, ignore.case = TRUE)
  
  # タイトルに除外キーワードが含まれ、かつ抄録に特定のキーワードが含まれない場合のみ除外
  title_excluded && !abstract_included
}

# Excl1列を追加し、除外基準に該当する論文にフラグを立てる
df <- df %>%
  mutate(Excl1 = ifelse(sapply(seq_len(nrow(df)), function(i) is_excluded(df$Title[i], df$Study[i])), 
                        paste(Author, Year, sep=" "),
                        NA))

# 更新したデータフレームをCSVファイルとして保存
write.csv(df, "dfMeta_updated.csv", row.names = FALSE)

# 除外された論文数と残った論文数を表示
excluded_count <- sum(!is.na(df$Excl1))
remaining_count <- nrow(df) - excluded_count

cat("除外された論文数:", excluded_count, "\n")
cat("残った論文数:", remaining_count, "\n")
### End of Escl1
#############################################
### Start of Excl2
# CSVファイルの読み込み
df <- read.csv("dfMeta_updated.csv", stringsAsFactors = FALSE)

# Excl2の除外基準（新しい基準を設定）
excl2_keywords <- c("animal", "in vitro", "cell culture", "molecular", "genetic")

# 除外基準の関数を定義
is_excluded <- function(title, keywords) {
  lower_title <- tolower(title)
  any(sapply(keywords, function(keyword) str_detect(lower_title, keyword)))
}

# Excl2列を追加（Excl1の内容を保持しつつ、新たな除外基準を適用）
df <- df %>%
  mutate(Excl2 = case_when(
    !is.na(Excl1) ~ NA_character_,  # Excl1に値がある場合はExcl2をNAに
    sapply(Title, function(t) is_excluded(t, excl2_keywords)) ~ paste(Author, Year, sep=" "),
    TRUE ~ NA_character_
  ))

# 重複文献の除外前の総数
total_before <- nrow(df)

# 重複文献の除外
df <- distinct(df, Title, .keep_all = TRUE)

# 重複除外による減少数
duplicates_removed <- total_before - nrow(df)

# 除外された論文数と残った論文数を表示
excluded_count_excl1 <- sum(!is.na(df$Excl1))
excluded_count_excl2 <- sum(!is.na(df$Excl2))
remaining_count <- nrow(df) - excluded_count_excl1 - excluded_count_excl2

cat("Excl1で除外された論文数:", excluded_count_excl1, "\n")
cat("Excl2で除外された論文数:", excluded_count_excl2, "\n")
cat("重複により除外された論文数:", duplicates_removed, "\n")
cat("残った論文数:", remaining_count, "\n")

# Excl2で除外されたタイトルを表示
# cat("\nExcl2で除外されたタイトル:\n")
# excluded_titles <- df$Title[!is.na(df$Excl2)]
# for (title in excluded_titles) {
#   cat(title, "\n")
# }

# 最終的な論文リストをCSVファイルとして保存
write.csv(df, "dfMeta_02.csv", row.names = FALSE)
### Excl2 End ###
#############################################
### Start of Excl3

library(rentrez)
library(dplyr)

# PMIDのリストを読み込む（dfMeta_cycling.csvから）
df_original <- read.csv("dfMeta_cycling.csv", stringsAsFactors = FALSE)

# Abstractを取得する関数（修正版）
get_abstract <- function(pmid) {
  tryCatch({
    summary <- entrez_summary(db="pubmed", id=pmid)
    abstract <- summary$abstract
    if (is.null(abstract)) {
      return(NA_character_)
    } else if (is.list(abstract)) {
      return(paste(unlist(abstract), collapse=" "))
    } else {
      return(as.character(abstract))
    }
  }, error = function(e) {
    return(NA_character_)
  })
}

# AbstractカラムをPMIDから取得して追加
df_updated <- df_original %>%
  mutate(Abstract = sapply(PMID, get_abstract))

# 更新したデータフレームをCSVファイルとして保存
write.csv(df_updated, "dfMeta_updated.csv", row.names = FALSE)
#############################################
library(dplyr)
library(stringr)

# CSVファイルの読み込み
df <- read.csv("dfMeta_updated.csv", stringsAsFactors = FALSE)

# ファイルの構造を確認
print(names(df))

# Excl1とExcl2カラムが存在しない場合は作成
if(!"Excl1" %in% names(df)) {
  df$Excl1 <- NA_character_
}
if(!"Excl2" %in% names(df)) {
  df$Excl2 <- NA_character_
}

# Excl3の除外基準を定義
is_excluded_excl3 <- function(title, abstract) {
  lower_title <- tolower(title)
  lower_abstract <- tolower(as.character(abstract))  # NAを処理するためにas.character()を使用
  
  # Abstractが NA の場合はタイトルのみを使用
  if (is.na(abstract)) {
    lower_abstract <- lower_title
  }
  
  # 1. 研究デザインが適切でない
  if (str_detect(lower_abstract, "cross-sectional|case report|case study")) {
    return("Inappropriate study design")
  }
  
  # 2. 対象集団が適切でない
  if (str_detect(lower_abstract, "patients|disease|disorder") & 
      !str_detect(lower_abstract, "healthy")) {
    return("Inappropriate population")
  }
  
  # 3. 介入が自転車運動に関連していない
  if (!str_detect(lower_abstract, "cycl|bik|pedal")) {
    return("Intervention not related to cycling")
  }
  
  # 4. アウトカムが健康・フィットネスに関連していない
  health_keywords <- c("health", "fitness", "cardiovascular", "metabolic", "strength", "endurance")
  if (!any(sapply(health_keywords, function(k) str_detect(lower_abstract, k)))) {
    return("Outcome not related to health/fitness")
  }
  
  return(NA_character_)
}

# Excl3列を追加
df <- df %>%
  mutate(Excl3 = case_when(
    !is.na(Excl1) | !is.na(Excl2) ~ NA_character_,  # Excl1またはExcl2に値がある場合はExcl3をNAに
    TRUE ~ sapply(1:n(), function(i) is_excluded_excl3(Title[i], Abstract[i]))
  ))

# 除外された論文数と残った論文数を表示
excluded_count_excl1 <- sum(!is.na(df$Excl1))
excluded_count_excl2 <- sum(!is.na(df$Excl2))
excluded_count_excl3 <- sum(!is.na(df$Excl3))
remaining_count <- nrow(df) - excluded_count_excl1 - excluded_count_excl2 - excluded_count_excl3

cat("Excl1で除外された論文数:", excluded_count_excl1, "\n")
cat("Excl2で除外された論文数:", excluded_count_excl2, "\n")
cat("Excl3で除外された論文数:", excluded_count_excl3, "\n")
cat("残った論文数:", remaining_count, "\n")

# Excl3で除外された理由を表示
cat("\nExcl3で除外された理由:\n")
excluded_reasons <- table(df$Excl3[!is.na(df$Excl3)])
print(excluded_reasons)

# 最終的な論文リストをCSVファイルとして保存
write.csv(df, "dfMeta_03.csv", row.names = FALSE)
### Excl3 End ###
#############################################
### Start of Excl4
library(dplyr)
library(stringr)

# CSVファイルの読み込み（前回の結果を使用）
df <- read.csv("dfMeta_03.csv", stringsAsFactors = FALSE)

# ROBINS-Iの各ドメインを評価する関数（より厳密な基準）
assess_robins_i <- function(title, abstract) {
  lower_text <- tolower(paste(title, abstract, sep = " "))
  
  score <- 0
  
  # 各ドメインのスコアリング
  if(str_detect(lower_text, "adjust|control|accounted for")) score <- score + 1
  if(str_detect(lower_text, "representative sample|population-based")) score <- score + 1
  if(str_detect(lower_text, "clear definition|well-defined intervention")) score <- score + 1
  if(str_detect(lower_text, "adherence|compliance")) score <- score + 1
  if(str_detect(lower_text, "complete data|no missing|imputation")) score <- score + 1
  if(str_detect(lower_text, "validated|objective measure")) score <- score + 1
  if(str_detect(lower_text, "pre-specified|registered protocol")) score <- score + 1
  
  # スコアに基づく評価
  if(score >= 5) {
    return("Low risk of bias")
  } else if(score >= 3) {
    return("Moderate risk of bias")
  } else {
    return("High risk of bias")
  }
}

# Excl4列を追加
df <- df %>%
  mutate(Excl4 = sapply(1:n(), function(i) assess_robins_i(Title[i], Abstract[i])))

# 高リスクと中リスクの研究を除外
df <- df %>%
  mutate(Excl4 = ifelse(Excl4 %in% c("High risk of bias", "Moderate risk of bias"), Excl4, NA_character_))

# 除外された論文数と残った論文数を表示
excluded_count_excl1 <- sum(!is.na(df$Excl1))
excluded_count_excl2 <- sum(!is.na(df$Excl2))
excluded_count_excl3 <- sum(!is.na(df$Excl3))
excluded_count_excl4 <- sum(!is.na(df$Excl4))
remaining_count <- nrow(df) - excluded_count_excl1 - excluded_count_excl2 - excluded_count_excl3 - excluded_count_excl4

cat("Excl1で除外された論文数:", excluded_count_excl1, "\n")
cat("Excl2で除外された論文数:", excluded_count_excl2, "\n")
cat("Excl3で除外された論文数:", excluded_count_excl3, "\n")
cat("Excl4で除外された論文数:", excluded_count_excl4, "\n")
cat("最終的に残った論文数:", remaining_count, "\n")

# Excl4での評価結果を表示
cat("\nROBINS-Iによる評価結果:\n")
robins_i_results <- table(df$Excl4, useNA = "ifany")
print(robins_i_results)

# 残った論文のタイトルを表示
cat("\n残った論文のタイトル:\n")
remaining_papers <- df %>% 
  filter(is.na(Excl1) & is.na(Excl2) & is.na(Excl3) & is.na(Excl4)) %>%
  select(Title)
print(remaining_papers)

# 最終的な論文リストをCSVファイルとして保存
write.csv(df, "dfMeta_04.csv", row.names = FALSE)

### Abstractの有無を確認する#############################################
# Abstractを取得する関数の改善
get_abstract <- function(pmid) {
  tryCatch({
    article <- entrez_fetch(db="pubmed", id=pmid, rettype="xml", parsed=TRUE)
    abstract <- xml_text(xml_find_first(article, ".//Abstract/AbstractText"))
    if(is.null(abstract) || abstract == "") {
      return(NA_character_)
    }
    return(abstract)
  }, error = function(e) {
    cat("Error fetching abstract for PMID", pmid, ":", conditionMessage(e), "\n")
    return(NA_character_)
  })
}

# データフレームの作成（Abstractの取得を含む）
df <- do.call(rbind, lapply(all_summaries, function(x) {
  abstract <- get_abstract(x$uid)
  data.frame(
    PMID = ifelse(is.null(x$uid), NA, x$uid),
    Author = ifelse(is.null(x$sortfirstauthor), NA, x$sortfirstauthor),
    Year = ifelse(is.null(x$pubdate), NA, substr(x$pubdate, 1, 4)),
    Title = ifelse(is.null(x$title), NA, x$title),
    Abstract = abstract,
    stringsAsFactors = FALSE
  )
}))

# CSVファイルとして保存（エンコーディングを指定）
write.csv(df, "dfMeta_final.csv", row.names = FALSE, fileEncoding = "UTF-8")

# 保存したファイルを読み込んで確認
df_check <- read.csv("dfMeta_final.csv", stringsAsFactors = FALSE, encoding = "UTF-8")

# Abstractの状態を確認
abstract_status <- table(is.na(df_check$Abstract))
print(abstract_status)
