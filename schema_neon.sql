--
-- PostgreSQL database dump
--

\restrict y5tGbYhgHv3p1c20v5exUb1HiYjYdvF1BmC2moMXdWu7VyHRdXDxoXEXQXTLtLQ

-- Dumped from database version 18.3 (Ubuntu 18.3-1.pgdg24.04+1)
-- Dumped by pg_dump version 18.3 (Ubuntu 18.3-1.pgdg24.04+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET transaction_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: cpg_annotation; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.cpg_annotation (
    site_id text NOT NULL,
    gene text,
    chr text,
    start_pos integer,
    end_pos integer,
    distance_tss text,
    cgi text,
    cg_island text
);


--
-- Name: cpg_features; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.cpg_features (
    site_id text NOT NULL,
    leukocyte_median real,
    leukocyte_std real,
    n_samples integer
);


--
-- Name: cpg_gene_map; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.cpg_gene_map (
    site_id text NOT NULL,
    ensembl_id text NOT NULL,
    gene_symbol text
);


--
-- Name: expression_correlation; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.expression_correlation (
    site_id text NOT NULL,
    tumor_type text NOT NULL,
    spearman_r real,
    n_samples integer
);


--
-- Name: gene_annotation; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.gene_annotation (
    ensembl_id text NOT NULL,
    gene_symbol text,
    chr text,
    start_pos integer,
    end_pos integer,
    strand text,
    biotype text,
    tss integer
);


--
-- Name: sample_metadata; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.sample_metadata (
    sample_id text NOT NULL,
    patient_id text,
    tumor_type text,
    sample_type text,
    sample_type_id integer,
    sex text,
    age integer,
    tissue_type text,
    sample_class text,
    platform text,
    batch text
);


--
-- Name: tumor_summary; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.tumor_summary (
    site_id text NOT NULL,
    tumor_type text NOT NULL,
    tumor_median real,
    tumor_std real,
    normal_median real,
    normal_std real,
    pan_tumor_median real,
    pan_tumor_std real,
    pan_normal_median real,
    pan_normal_std real,
    delta_median real,
    dispersion_index real,
    n_tumor integer,
    n_normal integer
);


--
-- Name: tumor_types; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.tumor_types (
    tumor_type text NOT NULL,
    full_name text,
    tissue text
);


--
-- Name: cpg_annotation cpg_annotation_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.cpg_annotation
    ADD CONSTRAINT cpg_annotation_pkey PRIMARY KEY (site_id);


--
-- Name: cpg_features cpg_features_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.cpg_features
    ADD CONSTRAINT cpg_features_pkey PRIMARY KEY (site_id);


--
-- Name: gene_annotation gene_annotation_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.gene_annotation
    ADD CONSTRAINT gene_annotation_pkey PRIMARY KEY (ensembl_id);


--
-- Name: sample_metadata sample_metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.sample_metadata
    ADD CONSTRAINT sample_metadata_pkey PRIMARY KEY (sample_id);


--
-- Name: tumor_summary tumor_summary_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tumor_summary
    ADD CONSTRAINT tumor_summary_pkey PRIMARY KEY (site_id, tumor_type);


--
-- Name: tumor_types tumor_types_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.tumor_types
    ADD CONSTRAINT tumor_types_pkey PRIMARY KEY (tumor_type);


--
-- Name: idx_cpg_annotation_chr_start_pos; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_cpg_annotation_chr_start_pos ON public.cpg_annotation USING btree (chr, start_pos);


--
-- Name: idx_cpg_ensembl; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_cpg_ensembl ON public.cpg_gene_map USING btree (ensembl_id);


--
-- Name: idx_cpg_gene_symbol; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_cpg_gene_symbol ON public.cpg_gene_map USING btree (gene_symbol);


--
-- Name: idx_expr_corr_rho; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_expr_corr_rho ON public.expression_correlation USING btree (spearman_r);


--
-- Name: idx_expr_corr_site; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_expr_corr_site ON public.expression_correlation USING btree (site_id);


--
-- Name: idx_expr_corr_tumor; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_expr_corr_tumor ON public.expression_correlation USING btree (tumor_type);


--
-- Name: idx_gene_annotation_symbol; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_gene_annotation_symbol ON public.gene_annotation USING btree (gene_symbol);


--
-- Name: idx_tumor_summary_delta; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_tumor_summary_delta ON public.tumor_summary USING btree (delta_median DESC);


--
-- Name: idx_tumor_summary_tumor; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX idx_tumor_summary_tumor ON public.tumor_summary USING btree (tumor_type);


--
-- PostgreSQL database dump complete
--

\unrestrict y5tGbYhgHv3p1c20v5exUb1HiYjYdvF1BmC2moMXdWu7VyHRdXDxoXEXQXTLtLQ

