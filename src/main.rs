pub mod fft;

use anyhow::{Result, anyhow};
use fft::fft;
use cpal::{
    traits::{DeviceTrait, HostTrait, StreamTrait}, Device, FromSample, Sample, SizedSample
};

fn get_default_device_config() -> Result<Device> {
    let host = cpal::default_host();
    let device = host.default_output_device();
    let device = device.ok_or(anyhow!("None device founded."))?;

    println!("Default host: {}", device.name().unwrap_or("null".to_string()));

    if let Ok(config) = device.default_output_config() {
        println!("SampleFormat: {}", config.sample_format());
        println!("SampleRate: {}", config.sample_rate().0);
        println!("Channels: {}", config.channels());
    }

    Ok(device)
}

fn main() {
    println!("Hello, world! {}...{}", u16::MIN, u16::MAX);

    let _ = get_default_device_config();
}
